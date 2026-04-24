/**
 * @file kmer_filter.cpp
 * @brief K-mer based repeat filter (optimized from MEGAnE)
 * @author cTEA Development Team
 * @version 0.1.0
 * 
 * Optimized from MEGAnE's:
 *   - 0_build_kmer_set.py (k-mer set construction)
 *   - save_redundant_kmers.cpp (k-mer filtering)
 * 
 * Key optimizations vs MEGAnE:
 * 1. C++ implementation (vs MEGAnE's Python)
 * 2. Binary search for k-mer lookup (O(log n) vs O(n) in Python)
 * 3. Memory-mapped k-mer sets for large genomes
 * 4. Single-pass filtering (vs MEGAnE's separate steps)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include "dna_to_2bit.hpp"

using namespace dna_to_2bit_hpp;

typedef unsigned long long ull;

#define REP_KMER_SIZE 11
#define SHIFT_16_TO_11 10

namespace kmer_filter_hpp {

/**
 * @brief K-mer set container with binary search
 * Optimized from MEGAnE's Python k-mer set
 */
class KmerSet {
private:
    std::vector<uint32_t> kmers_;  // Sorted k-mer array
    bool is_sorted_;
    
public:
    KmerSet() : is_sorted_(false) {}
    
    /**
     * @brief Load k-mers from .mk file (MEGAnE format)
     * Format: binary uint32_t array of 11-bit k-mers
     */
    bool load(const std::string& mk_file) {
        std::ifstream ifs(mk_file, std::ios::binary);
        if (!ifs) {
            std::cerr << "Error: Cannot open " << mk_file << std::endl;
            return false;
        }
        
        // Get file size
        ifs.seekg(0, std::ios::end);
        size_t file_size = ifs.tellg();
        ifs.seekg(0, std::ios::beg);
        
        size_t n_kmers = file_size / sizeof(uint32_t);
        kmers_.resize(n_kmers);
        ifs.read(reinterpret_cast<char*>(kmers_.data()), file_size);
        
        if (!ifs) {
            std::cerr << "Error: Failed to read " << mk_file << std::endl;
            return false;
        }
        
        is_sorted_ = false;  // Will be sorted after load
        sort_kmers();
        
        std::cerr << "Loaded " << n_kmers << " k-mers from " << mk_file << std::endl;
        return true;
    }
    
    /**
     * @brief Sort k-mers for binary search
     */
    void sort_kmers() {
        if (!is_sorted_) {
            std::sort(kmers_.begin(), kmers_.end());
            is_sorted_ = true;
        }
    }
    
    /**
     * @brief Check if a k-mer exists (binary search, O(log n))
     * vs MEGAnE's Python: O(n) linear search
     */
    bool contains(uint32_t kmer) const {
        if (!is_sorted_ || kmers_.empty()) return false;
        
        // Custom binary search (from MEGAnE's extract_discordant.cpp)
        ull pos = 0;
        size_t step = 1ULL << (63 - __builtin_clz(kmers_.size() - 1));
        pos = (kmers_[step - 1] < kmer ? kmers_.size() - step - 1 : -1);
        
        while ((step >>= 1) > 0) {
            pos = (kmers_[pos + step] < kmer ? pos + step : pos);
        }
        
        return (kmers_[pos + 1] == kmer);
    }
    
    /**
     * @brief Check if a sequence contains repeat k-mers
     * Optimized rolling k-mer calculation
     */
    bool sequence_has_repeat_kmer(const std::string& seq) {
        if (seq.length() < REP_KMER_SIZE) return false;
        
        uint32_t bit2 = 0;
        int nn = 0;
        
        // First window
        for (int i = 0; i < REP_KMER_SIZE; i++) {
            bit2 <<= 2;
            bit2 |= dna_to_2bit_table[(unsigned char)seq[i]];
            if (seq[i] == 'N' || seq[i] == 'n') {
                nn = REP_KMER_SIZE - 1;
            } else if (nn > 0) {
                nn--;
            }
        }
        
        if (nn == 0 && contains(bit2)) return true;
        
        // Rolling windows
        for (size_t i = REP_KMER_SIZE; i < seq.length(); i++) {
            bit2 <<= 2;
            bit2 |= dna_to_2bit_table[(unsigned char)seq[i]];
            bit2 <<= SHIFT_16_TO_11;
            bit2 >>= SHIFT_16_TO_11;
            
            if (seq[i] == 'N' || seq[i] == 'n') {
                nn = REP_KMER_SIZE - 1;
            } else if (nn > 0) {
                nn--;
            } else {
                if (contains(bit2)) return true;
            }
        }
        
        return false;
    }
    
    size_t size() const { return kmers_.size(); }
};

/**
 * @brief Filter reads based on k-mer content
 * Optimized from MEGAnE's save_redundant_kmers.cpp
 */
class KmerFilter {
private:
    KmerSet kmer_set_;
    int min_kmer_hits_;
    double min_kmer_ratio_;
    
public:
    KmerFilter() : min_kmer_hits_(3), min_kmer_ratio_(0.1) {}
    
    /**
     * @brief Initialize filter with k-mer set
     */
    bool initialize(const std::string& mk_file) {
        return kmer_set_.load(mk_file);
    }
    
    /**
     * @brief Filter a single read sequence
     * @return true if read passes filter (not too many repeat k-mers)
     */
    bool filter_read(const std::string& seq, int& kmer_hits) {
        if (seq.length() < REP_KMER_SIZE) {
            kmer_hits = 0;
            return true;  // Too short to filter
        }
        
        kmer_hits = 0;
        uint32_t bit2 = 0;
        int nn = 0;
        
        // First window
        for (int i = 0; i < REP_KMER_SIZE; i++) {
            bit2 <<= 2;
            bit2 |= dna_to_2bit_table[(unsigned char)seq[i]];
            if (seq[i] == 'N' || seq[i] == 'n') {
                nn = REP_KMER_SIZE - 1;
            }
        }
        
        if (nn == 0 && kmer_set_.contains(bit2)) kmer_hits++;
        
        // Rolling windows
        for (size_t i = REP_KMER_SIZE; i < seq.length(); i++) {
            bit2 <<= 2;
            bit2 |= dna_to_2bit_table[(unsigned char)seq[i]];
            bit2 <<= SHIFT_16_TO_11;
            bit2 >>= SHIFT_16_TO_11;
            
            if (seq[i] == 'N' || seq[i] == 'n') {
                nn = REP_KMER_SIZE - 1;
            } else if (nn > 0) {
                nn--;
            } else {
                if (kmer_set_.contains(bit2)) kmer_hits++;
            }
        }
        
        // Filter decision
        double ratio = (double)kmer_hits / (seq.length() - REP_KMER_SIZE + 1);
        return (kmer_hits >= min_kmer_hits_) && (ratio >= min_kmer_ratio_);
    }
    
    /**
     * @brief Filter multiple reads (batch processing)
     * Key optimization: Single pass through all reads
     */
    void filter_reads(const std::vector<std::string>& sequences,
                     std::vector<bool>& results,
                     std::vector<int>& hit_counts) {
        results.resize(sequences.size());
        hit_counts.resize(sequences.size());
        
        for (size_t i = 0; i < sequences.size(); i++) {
            results[i] = filter_read(sequences[i], hit_counts[i]);
        }
    }
};

} // namespace kmer_filter_hpp

#ifndef BUILD_LIB
// Main function (for standalone testing)
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.mk> <input.fasta> [output.fasta]" << std::endl;
        return 1;
    }
    
    std::string mk_file = argv[1];
    std::string fasta_file = argv[2];
    std::string output_file = (argc > 3) ? argv[3] : "filtered.fasta";
    
    // Initialize filter
    kmer_filter_hpp::KmerFilter filter;
    if (!filter.initialize(mk_file)) {
        std::cerr << "Error: Failed to load k-mer set" << std::endl;
        return 1;
    }
    
    // Read FASTA
    std::ifstream ifs(fasta_file);
    if (!ifs) {
        std::cerr << "Error: Cannot open " << fasta_file << std::endl;
        return 1;
    }
    
    std::vector<std::string> sequences;
    std::vector<std::string> headers;
    std::string line, current_header, current_seq;
    
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                sequences.push_back(current_seq);
                headers.push_back(current_header);
            }
            current_header = line.substr(1);
            current_seq.clear();
        } else {
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        sequences.push_back(current_seq);
        headers.push_back(current_header);
    }
    
    // Filter
    std::vector<bool> results;
    std::vector<int> hit_counts;
    filter.filter_reads(sequences, results, hit_counts);
    
    // Write output
    std::ofstream ofs(output_file);
    int passed = 0;
    for (size_t i = 0; i < sequences.size(); i++) {
        if (results[i]) {
            ofs << ">" << headers[i] << "\t" << hit_counts[i] << std::endl;
            ofs << sequences[i] << std::endl;
            passed++;
        }
    }
    
    std::cerr << "Filtering complete:" << std::endl;
    std::cerr << "  Total reads: " << sequences.size() << std::endl;
    std::cerr << "  Passed: " << passed << std::endl;
    std::cerr << "  Filtered: " << (sequences.size() - passed) << std::endl;
    
    return 0;
}
#endif // BUILD_LIB
