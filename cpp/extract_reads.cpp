/**
 * @file extract_reads.cpp
 * @brief Combined extractor for discordant and unmapped reads (optimized from MEGAnE)
 * @author cTEA Development Team
 * @version 0.2.0
 * 
 * Optimized from MEGAnE's:
 *   - extract_discordant.cpp (discordant pair detection)
 *   - extract_unmapped.cpp (unmapped read extraction)
 * 
 * Key optimizations vs MEGAnE:
 * 1. Single-pass scanning (vs MEGAnE's separate scans)
 * 2. In-memory data structures (vs MEGAnE's file-per-thread output)
 * 3. Avoids xTEA's chromosome-splitting bottleneck
 * 4. Direct output to Python-accessible format (for pybind11/ctypes)
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <functional>
#include <unistd.h>
#include <htslib/sam.h>
#include "dna_to_2bit.hpp"
#include "complementary_seq.hpp"
#include "bam_processor.h"
#include "ThreadPool.h"

using namespace dna_to_2bit_hpp;
using namespace complementary_seq_hpp;

#define MAX_CIGAR_LEN 128
#define TMP_BUF_SIZE  131072
#define READ_PAIR_GAP_LEN 2000
#define MAX_SEQ_LEN 512
#define DISCORDANT_READS_CLIP_LEN 20
#define UNMAPPED_MIN_LEN 13
#define REP_KMER_SIZE 11
#define SHIFT_16_TO_11 10
#define ABS_MIN_DIST 50
#define ABS_MAX_DIST 20000
#define MAPPED_REGION_LOW_COMPLEX_THRESHOLD 0.7
#define POLYA_OVERHANG_THRESHOLD 0.7

typedef unsigned long ul;
typedef unsigned long long ull;

const char ATGC[] = "ATGC";
const std::string SEMICOLON_STR = ";";

namespace extract_reads_hpp {

/**
 * @brief Structure to hold soft clip information from SA tags
 */
struct SoftClipInfo {
    std::string chr;
    uint32_t pos;
    char breakpoint;  // 'L' or 'R'
    int32_t clipstart;
    int32_t clipend;
    int32_t l_clip_len;
    int32_t r_clip_len;
    bool is_reverse;
};

/**
 * @brief Structure to hold extraction results in memory
 * Avoids MEGAnE's file-per-thread approach
 */
struct ExtractionResults {
    // pA reads: "chr:pos-end/strand/readname/is_read2"
    std::vector<std::string> pA_reads;
    
    // Overhang reads: "chr:pos-end/breakpoint/readname/is_read2"
    std::vector<std::string> overhang_reads;
    
    // Mapped regions: "seq/header" pairs
    std::unordered_map<std::string, std::string> mapped_regions;
    
    // Distant (hybrid) reads: "chr:pos-end/breakpoint/readname/is_read2"
    std::vector<std::string> distant_reads;
    
    // Absent reads: "readname/is_read2\tchr\tpos\tchr:pos-end\tstrand..."
    std::vector<std::string> absent_reads;
    
    // Unmapped reads: ">readname/is_read2\nseq"
    std::vector<std::string> unmapped_reads;
    
    // Statistics
    int64_t pA_count = 0;
    int64_t chimeric_count = 0;
    int64_t distant_count = 0;
    int64_t absent_count = 0;
    int64_t unmapped_count = 0;
    int64_t read_count = 0;
    
    void clear() {
        pA_reads.clear();
        overhang_reads.clear();
        mapped_regions.clear();
        distant_reads.clear();
        absent_reads.clear();
        unmapped_reads.clear();
        pA_count = chimeric_count = distant_count = absent_count = unmapped_count = read_count = 0;
    }
};

/**
 * @brief Main extractor class (optimized single-pass approach)
 */
class ReadExtractor {
private:
    samFile* bam_file_;
    bam_hdr_t* bam_header_;
    const std::vector<uint32_t>* crepkmer_;
    ull num_kmer_;
    const std::unordered_map<std::string, bool>* is_mainchr_;
    int n_threads_;
    bool is_cram_;
    std::string ref_fa_;
    
    // Results storage (in-memory, vs MEGAnE's file output)
    ExtractionResults results_;
    
    // Thread pool for parallel processing
    ThreadPool* thread_pool_;
    
public:
    ReadExtractor() : bam_file_(nullptr), bam_header_(nullptr), 
                     crepkmer_(nullptr), num_kmer_(0), 
                     is_mainchr_(nullptr), n_threads_(1), 
                     is_cram_(false), thread_pool_(nullptr) {}
    
    ~ReadExtractor() {
        if (thread_pool_) delete thread_pool_;
        if (bam_header_) bam_hdr_destroy(bam_header_);
        if (bam_file_) sam_close(bam_file_);
    }
    
    /**
     * @brief Initialize extractor with input file
     * Key optimization: Single file handle for all extractions
     */
    bool initialize(const std::string& bam_path, 
                   const std::string& ref_fa,
                   const std::vector<uint32_t>& rep_kmer,
                   ull num_kmer,
                   const std::unordered_map<std::string, bool>& mainchrs,
                   int n_threads) {
        // Check if CRAM
        is_cram_ = (bam_path.find(".cram") != std::string::npos);
        ref_fa_ = ref_fa;
        
        // Open file
        bam_file_ = sam_open(bam_path.c_str(), "r");
        
        // For CRAM, set reference if provided
        if (is_cram_ && !ref_fa.empty()) {
            if (hts_set_fai_filename(bam_file_, ref_fa.c_str()) != 0) {
                std::cerr << "Warning: Failed to set reference for CRAM" << std::endl;
            }
        }
        
        if (!bam_file_) {
            std::cerr << "Error: Cannot open " << bam_path << std::endl;
            return false;
        }
        
        bam_header_ = sam_hdr_read(bam_file_);
        if (!bam_header_) {
            std::cerr << "Error: Cannot read header" << std::endl;
            return false;
        }
        
        crepkmer_ = &rep_kmer;
        num_kmer_ = num_kmer;
        is_mainchr_ = &mainchrs;
        n_threads_ = n_threads;
        
        if (n_threads_ > 1) {
            thread_pool_ = new ThreadPool(n_threads_);
        }
        
        return true;
    }
    
    /**
     * @brief Single-pass extraction (KEY OPTIMIZATION vs MEGAnE)
     */
    bool extract_all(ExtractionResults& results) {
        results_ = results;  // Reference to output
        
        bam1_t* aln = bam_init1();
        int ret;
        int64_t processed = 0;
        
        while ((ret = sam_read1(bam_file_, bam_header_, aln)) >= 0) {
            process_alignment(aln);
            bam_destroy1(aln);
            aln = bam_init1();
            
            processed++;
            if (processed % 1000000 == 0) {
                std::cerr << "Processed " << processed << " reads..." << std::endl;
            }
        }
        
        if (aln) bam_destroy1(aln);
        
        results = results_;  // Copy results back
        return true;
    }
    
    /**
     * @brief Process a single alignment
     */
    void process_alignment(bam1_t* aln) {
        results_.read_count++;
        
        uint16_t flag = aln->core.flag;
        
        // Skip supplementary alignments
        if (flag & BAM_FSUPPLEMENTARY) return;
        
        // Skip single-end reads
        if (!(flag & BAM_FPAIRED)) return;
        
        // === Case 1: Unmapped read ===
        if (flag & BAM_FUNMAP) {
            extract_unmapped(aln);
            return;
        }
        
        // === Case 2: Discordant pair or chimeric read ===
        int32_t isize = aln->core.isize;
        bool is_discordant = (isize == 0) || (isize <= -READ_PAIR_GAP_LEN) || (isize >= READ_PAIR_GAP_LEN);
        
        // Check for SA tag
        uint8_t* sa_tag = bam_aux_get(aln, "SA");
        bool has_sa = (sa_tag != nullptr);
        
        if (is_discordant || has_sa) {
            extract_discordant(aln, has_sa, sa_tag);
        }
    }
    
    /**
     * @brief Extract unmapped reads
     */
    void extract_unmapped(bam1_t* aln) {
        char* qname = bam_get_qname(aln);
        int32_t l_qseq = aln->core.l_qseq;
        uint8_t* seq_p = bam_get_seq(aln);
        
        if (l_qseq < UNMAPPED_MIN_LEN) return;
        
        std::string seq;
        seq.reserve(l_qseq);
        int non_N_count = 0;
        for (int i = 0; i < l_qseq; i++) {
            char nt = seq_nt16_str[bam_seqi(seq_p, i)];
            if (nt != 'N' && nt != 'n') {
                non_N_count++;
            }
            seq += nt;
        }
        
        if (non_N_count < UNMAPPED_MIN_LEN) return;
        
        if (is_rep_kmer(seq, l_qseq)) {
            bool is_read2 = (aln->core.flag & BAM_FREAD2) > 0;
            std::string header = std::string(qname) + "/" + std::to_string((int)is_read2 + 1);
            results_.unmapped_reads.push_back(">" + header + "\n" + seq + "\n");
            results_.unmapped_count++;
        }
    }
    
    /**
     * @brief Extract discordant/chimeric reads
     */
    void extract_discordant(bam1_t* aln, bool has_sa, uint8_t* sa_tag) {
        uint32_t* cigar = bam_get_cigar(aln);
        uint32_t n_cigar = aln->core.n_cigar;
        
        bool has_H = false, has_S = false;
        std::pair<char, uint32_t> cigar_arr[MAX_CIGAR_LEN];
        parse_cigar(cigar, n_cigar, cigar_arr, has_H, has_S);
        
        if (has_H) return;
        
        if (!has_S && !has_sa) return;
        
        int32_t l_clip_len = 0, r_clip_len = 0;
        int64_t clipstart = 0, clipend = 0;
        char breakpoint = 'N';
        int32_t l_qseq = aln->core.l_qseq;
        define_breakpoint(cigar_arr, n_cigar, breakpoint, l_clip_len, r_clip_len, clipstart, clipend);
        
        if (breakpoint == 'N') return;
        
        if (l_clip_len < DISCORDANT_READS_CLIP_LEN && r_clip_len < DISCORDANT_READS_CLIP_LEN) {
            return;
        }
        
        uint8_t* seq_p = bam_get_seq(aln);
        std::string fseq, rseq;
        fseq.reserve(l_qseq);
        rseq.reserve(l_qseq);
        for (int i = 0; i < l_qseq; i++) {
            char nt = seq_nt16_str[bam_seqi(seq_p, i)];
            fseq += nt;
            rseq += complementary_seq_hpp::complement(nt);
        }
        std::reverse(rseq.begin(), rseq.end());
        
        bool is_reverse = (aln->core.flag & BAM_FREVERSE) > 0;
        std::string chr = bam_header_->target_name[aln->core.tid];
        
        if (!is_mainchr_->at(chr)) return;
        
        if (has_sa) {
            std::vector<SoftClipInfo> softclips;
            char* qname = bam_get_qname(aln);
            parse_SA_tag((char*)(sa_tag + 1), softclips, cigar_arr, cigar, n_cigar,
                         breakpoint, l_clip_len, r_clip_len, clipstart, clipend, l_qseq);
            
            check_pA_or_overhang(softclips, fseq, rseq, l_qseq, is_reverse, chr, 
                               qname, l_clip_len, r_clip_len, clipstart, clipend);
        }
        
        char* qname = bam_get_qname(aln);
        check_overhang(fseq, rseq, l_qseq, is_reverse, chr, qname, 
                       l_clip_len, r_clip_len, clipstart, clipend);
    }
    
    void check_pA_or_overhang(const std::vector<SoftClipInfo>& softclips,
                              const std::string& fseq, const std::string& rseq,
                              int32_t l_qseq, bool is_reverse,
                              const std::string& chr, const char* qname,
                              int32_t l_clip_len, int32_t r_clip_len,
                              int64_t clipstart, int64_t clipend) {
        for (const auto& sc : softclips) {
            if (sc.breakpoint == 'N') continue;
            if (!is_mainchr_->at(sc.chr)) continue;
            
            std::string seq_to_use;
            if (sc.is_reverse && is_reverse) {
                seq_to_use = fseq;
            } else if (!sc.is_reverse && !is_reverse) {
                seq_to_use = fseq;
            } else {
                seq_to_use = rseq;
            }
            
            std::string clipped_seq;
            if (sc.breakpoint == 'L') {
                clipped_seq = seq_to_use.substr(sc.clipstart, sc.clipend - sc.clipstart);
            } else {
                clipped_seq = seq_to_use.substr(sc.clipstart);
            }
            
            int64_t Acount = 0;
            for (char c : clipped_seq) {
                if (c == 'A' || c == 'T') Acount++;
            }
            
            bool is_pA = ((double)Acount / clipped_seq.length()) >= POLYA_OVERHANG_THRESHOLD;
            
            if (is_pA) {
                char strand = sc.is_reverse ? '-' : '+';
                std::string header = std::string(qname) + "/" + 
                                    std::to_string((int)(sc.breakpoint == 'L' ? 1 : 2)) + 
                                    "/" + strand;
                results_.pA_reads.push_back(sc.chr + ":" + std::to_string(sc.pos) + "/" + 
                                           std::to_string(sc.clipend - sc.clipstart) + "/" + 
                                           header);
                results_.pA_count++;
            } else {
                if (is_rep_kmer(clipped_seq, clipped_seq.length())) {
                    char strand = sc.is_reverse ? '-' : '+';
                    std::string header = std::string(qname) + "/" + 
                                        std::to_string((int)(sc.breakpoint == 'L' ? 1 : 2)) + 
                                        "/" + strand;
                    results_.overhang_reads.push_back(sc.chr + ":" + std::to_string(sc.pos) + 
                                                      "/" + std::to_string(sc.breakpoint) + "/" + 
                                                      header);
                    results_.chimeric_count++;
                }
            }
            
            std::string mapped_seq = seq_to_use.substr(sc.l_clip_len, 
                                                        l_qseq - sc.l_clip_len - sc.r_clip_len);
            if (!is_simple_repeat(mapped_seq, mapped_seq.length())) {
                std::string mapped_header = std::string(qname) + "/" + 
                                            std::to_string((int)(sc.breakpoint == 'L' ? 1 : 2));
                results_.mapped_regions[mapped_seq] += mapped_header + ";";
            }
        }
    }
    
    void check_overhang(const std::string& fseq, const std::string& rseq,
                       int32_t l_qseq, bool is_reverse,
                       const std::string& chr, const char* qname,
                       int32_t l_clip_len, int32_t r_clip_len,
                       int64_t clipstart, int64_t clipend) {
        // Simplified for brevity
    }
    
    inline void parse_cigar(uint32_t* cigar, uint32_t& n_cigar,
                            std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                            bool& contains_H, bool& contains_S) {
        for (int i = 0; i < n_cigar; ++i) {
            char opchr = bam_cigar_opchr(cigar[i]);
            cigar_arr[i] = std::make_pair(opchr, bam_cigar_oplen(cigar[i]));
            if (opchr == 'H') contains_H = true;
            if (opchr == 'S') contains_S = true;
        }
    }
    
    inline void define_breakpoint(std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN], 
                                  uint32_t& n_cigar, char& breakpoint, 
                                  int32_t& l_clip_len, int32_t& r_clip_len,
                                  int64_t& clipstart, int64_t& clipend) {
        if (cigar_arr[0].first == 'S' || cigar_arr[0].first == 'H') {
            l_clip_len = cigar_arr[0].second;
        } else {
            l_clip_len = 0;
        }
        
        if (cigar_arr[n_cigar - 1].first == 'S' || cigar_arr[n_cigar - 1].first == 'H') {
            r_clip_len = cigar_arr[n_cigar - 1].second;
        } else {
            r_clip_len = 0;
        }
        
        if (l_clip_len == r_clip_len) {
            breakpoint = 'N';
            clipstart = 0;
            clipend = 0;
        } else if (l_clip_len > r_clip_len) {
            breakpoint = 'L';
            clipstart = 0;
            clipend = l_clip_len;
        } else {
            breakpoint = 'R';
            clipstart = 0;
            clipend = 0;
        }
    }
    
    inline bool is_simple_repeat(const std::string& seq, int32_t seqlen) {
        int32_t count = 0;
        for (char c : ATGC) {
            for (int32_t i = 0; i < seqlen; i++) {
                if (seq[i] == c) count++;
            }
            if (((double)count / seqlen) >= MAPPED_REGION_LOW_COMPLEX_THRESHOLD) {
                return true;
            }
            count = 0;
        }
        return false;
    }
    
    inline bool is_rep_kmer(const std::string& seq, uint64_t clip_len) {
        return false;  // Placeholder
    }
    
    void parse_SA_tag(char* sa_tag, std::vector<SoftClipInfo>& softclips,
                      std::pair<char, uint32_t> (&cigar_arr)[MAX_CIGAR_LEN],
                      uint32_t* cigar, uint32_t& n_cigar,
                      char& breakpoint, int32_t& l_clip_len, int32_t& r_clip_len,
                      int64_t& clipstart, int64_t& clipend, int32_t& l_qseq) {
        // Placeholder for SA tag parsing
    }
};

} // namespace extract_reads_hpp

#ifndef BUILD_LIB
// Main function (for standalone testing)
int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <input.bam/cram> <main_chrs.txt> <input.mk> <output_dir> [n_thread] [reference.fa]" << std::endl;
        return 1;
    }
    
    std::string bam_path = argv[1];
    std::string mainchr_file = argv[2];
    std::string mk_file = argv[3];
    std::string outdir = argv[4];
    int n_thread = (argc > 5) ? atoi(argv[5]) : 1;
    std::string ref_fa = (argc > 6) ? argv[6] : "";
    
    std::unordered_map<std::string, bool> mainchrs;
    std::vector<uint32_t> rep_kmer;
    ull num_kmer = 0;
    
    extract_reads_hpp::ReadExtractor extractor;
    if (!extractor.initialize(bam_path, ref_fa, rep_kmer, num_kmer, mainchrs, n_thread)) {
        std::cerr << "Error: Initialization failed" << std::endl;
        return 1;
    }
    
    extract_reads_hpp::ExtractionResults results;
    if (!extractor.extract_all(results)) {
        std::cerr << "Error: Extraction failed" << std::endl;
        return 1;
    }
    
    std::cerr << "Extraction complete:" << std::endl;
    std::cerr << "  pA reads: " << results.pA_count << std::endl;
    std::cerr << "  Chimeric reads: " << results.chimeric_count << std::endl;
    std::cerr << "  Distant reads: " << results.distant_count << std::endl;
    std::cerr << "  Absent reads: " << results.absent_count << std::endl;
    std::cerr << "  Unmapped reads: " << results.unmapped_count << std::endl;
    std::cerr << "  Total reads processed: " << results.read_count << std::endl;
    
    return 0;
}
#endif // BUILD_LIB