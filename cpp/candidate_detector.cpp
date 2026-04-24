/**
 * @file candidate_detector.cpp
 * @brief Candidate MEI detector (MAJOR OPTIMIZATION vs xTEA)
 * @author cTEA Development Team
 * @version 0.1.0
 * 
 * KEY OPTIMIZATION vs xTEA:
 *   xTEA: for chrm in chrms: for each BAM: read temp files, filter by chrm
 *         (x_TEI_locator.py lines 112-168)
 *   cTEA: SINGLE PASS -> populate hash maps in memory -> output
 *         (no intermediate files, no chromosome splitting)
 * 
 * Combines:
 *   - xTEA's clip+disc logic (from x_TEI_locator.py)
 *   - MEGAnE's speed (C++ implementation)
 *   - In-memory evidence accumulation (vs file-based)
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <functional>
#include <cmath>
#include <htslib/sam.h>
#include "bam_processor.h"

namespace candidate_detector_hpp {

/**
 * @brief Candidate site with accumulated evidence
 * Replaces xTEA's chromosome-splitting approach
 */
struct CandidateSite {
    std::string chrom;
    uint32_t pos;              // 0-based position
    
    // Evidence counts
    uint32_t left_clip_count;
    uint32_t right_clip_count;
    uint32_t left_clip_consensus;  // clipped part aligned to repeat consensus
    uint32_t right_clip_consensus;
    
    uint32_t left_discordant;
    uint32_t right_discordant;
    
    // PolyA evidence
    uint32_t left_polyA;
    uint32_t right_polyA;
    
    // Coverage info
    float left_coverage;
    float right_coverage;
    
    // Additional info
    std::string te_type;          // L1, Alu, SVA, etc.
    std::string subfamily;        // Specific subfamily
    float divergence;              // RepeatMasker divergence
    
    // Consensus alignment positions
    uint32_t left_cns_start;
    uint32_t left_cns_end;
    uint32_t right_cns_start;
    uint32_t right_cns_end;
    
    // Cluster info
    uint32_t left_cluster_start;
    uint32_t left_cluster_end;
    uint32_t right_cluster_start;
    uint32_t right_cluster_end;
    
    // Support type (from xTEA's classification)
    std::string support_type;  // two_side, one_half, one_side, other
    
    CandidateSite() : chrom(""), pos(0),
        left_clip_count(0), right_clip_count(0),
        left_clip_consensus(0), right_clip_consensus(0),
        left_discordant(0), right_discordant(0),
        left_polyA(0), right_polyA(0),
        left_coverage(0.0), right_coverage(0.0),
        divergence(100.0), left_cns_start(0), left_cns_end(0),
        right_cns_start(0), right_cns_end(0),
        left_cluster_start(0), left_cluster_end(0),
        right_cluster_start(0), right_cluster_end(0) {}
};

/**
 * @brief Main candidate detector class
 * SINGLE-PASS vs xTEA's chromosome-splitting
 */
class CandidateDetector {
private:
    // In-memory storage (vs xTEA's temp files)
    // Key: "chr:pos", Value: CandidateSite
    std::unordered_map<std::string, CandidateSite> candidates_;
    
    // BAM processor
    cTEA::BamProcessor* bam_processor_;
    
    // Parameters (from xTEA's global_values.py)
    uint32_t min_clip_len_;
    uint32_t min_disc_count_;
    float min_af_cutoff_;
    uint32_t clip_cluster_diff_cutoff_;
    uint32_t sva_clip_cluster_diff_cutoff_;
    
    // Statistics
    uint64_t total_reads_processed_;
    uint64_t candidates_found_;
    
public:
    CandidateDetector() : bam_processor_(nullptr),
        min_clip_len_(20),
        min_disc_count_(3),
        min_af_cutoff_(0.075),
        clip_cluster_diff_cutoff_(300),
        sva_clip_cluster_diff_cutoff_(200),
        total_reads_processed_(0),
        candidates_found_(0) {}
    
    ~CandidateDetector() {
        if (bam_processor_) {
            delete bam_processor_;
        }
    }
    
    /**
     * @brief Initialize detector
     */
    bool initialize(const std::string& bam_path,
                      const std::string& ref_path,
                      int n_threads) {
        bam_processor_ = new cTEA::BamProcessor();
        return bam_processor_->initialize(bam_path, ref_path, n_threads);
    }
    
    /**
     * @brief Main detection function (SINGLE-PASS)
     * 
     * This is the KEY OPTIMIZATION vs xTEA:
     * xTEA (x_TEI_locator.py lines 112-168):
     *   for chrm in m_chrms:
     *       for i in range(cnt):  # each BAM file
     *           if tmp_chrm != chrm: continue  # Skip if not matching
     * 
     * cTEA:
     *   Single scan -> populate hash maps -> output
     *   (No intermediate files, no chromosome iteration)
     */
    bool detect_candidates(
        const std::string& output_bed,
        const std::string& rep_kmer_file = "",
        const std::string& repeat_masker_bed = "") {
        
        if (!bam_processor_) {
            std::cerr << "Error: Detector not initialized" << std::endl;
            return false;
        }
        
        std::cerr << "Starting single-pass candidate detection..." << std::endl;
        std::cerr << "This replaces xTEA's chromosome-splitting approach" << std::endl;
        
        // Scan BAM and collect evidence in memory
        std::unordered_map<std::string, 
                          std::unordered_map<uint32_t, cTEA::CandidateEvidence>> evidence;
        
        auto progress_callback = [](uint32_t current, uint32_t total) {
            if (total == 0) {
                std::cerr << "\rProcessed " << current << " reads..." << std::flush;
            } else {
                std::cerr << "\rProcessed " << current << "/" << total 
                          << " reads (" << (current * 100 / total) << "%)..." << std::flush;
            }
        };
        
        if (!bam_processor_->scan_all_evidence(evidence, progress_callback)) {
            std::cerr << "\nError: BAM scan failed" << std::endl;
            return false;
        }
        
        std::cerr << "\nBAM scan complete. Building candidate sites..." << std::endl;
        
        // Convert evidence to candidate sites
        build_candidates(evidence);
        
        std::cerr << "Found " << candidates_.size() << " candidate sites" << std::endl;
        
        // Apply filters (merged from xTEA's logic)
        filter_candidates();
        
        // Output to BED file
        if (!output_bed.empty()) {
            return write_output_bed(output_bed);
        }
        
        return true;
    }
    
    /**
     * @brief Build candidate sites from evidence
     * Mimics xTEA's clustering logic but in memory
     */
    void build_candidates(
        const std::unordered_map<std::string, 
                      std::unordered_map<uint32_t, cTEA::CandidateEvidence>>& evidence) {
        
        for (const auto& chr_entry : evidence) {
            const std::string& chrom = chr_entry.first;
            const auto& pos_map = chr_entry.second;
            
            for (const auto& pos_entry : pos_map) {
                uint32_t pos = pos_entry.first;
                const auto& ev = pos_entry.second;
                
                // Create candidate
                CandidateSite site;
                site.chrom = chrom;
                site.pos = pos;
                
                // Copy evidence
                site.left_clip_count = ev.left_clip_count;
                site.right_clip_count = ev.right_clip_count;
                site.left_clip_consensus = ev.left_clip_consensus;
                site.right_clip_consensus = ev.right_clip_consensus;
                site.left_discordant = ev.left_discordant;
                site.right_discordant = ev.right_discordant;
                site.left_polyA = ev.left_polyA;
                site.right_polyA = ev.right_polyA;
                site.left_coverage = ev.left_coverage;
                site.right_coverage = ev.right_coverage;
                
                // Determine support type (from xTEA's logic)
                determine_support_type(site);
                
                // Only keep sites with enough evidence
                if (site.left_clip_count + site.right_clip_count >= min_clip_len_ ||
                    site.left_discordant + site.right_discordant >= min_disc_count_) {
                    std::string key = chrom + ":" + std::to_string(pos);
                    candidates_[key] = site;
                }
            }
        }
    }
    
    /**
     * @brief Determine support type (from xTEA's classification)
     * Mimics xTEA's get_ins_sub_type logic
     */
    void determine_support_type(CandidateSite& site) {
        bool has_left = (site.left_clip_count > 0 || site.left_discordant > 0);
        bool has_right = (site.right_clip_count > 0 || site.right_discordant > 0);
        
        if (has_left && has_right) {
            site.support_type = "two_side";
        } else if (has_left || has_right) {
            // Check if it's one_half (one side + half evidence on other side)
            if ((site.left_clip_count > 0 && site.right_discordant > 0) ||
                (site.right_clip_count > 0 && site.left_discordant > 0)) {
                site.support_type = "one_half";
            } else {
                site.support_type = "one_side";
            }
        } else {
            site.support_type = "other";
        }
    }
    
    /**
     * @brief Filter candidates (merged xTEA + MEGAnE logic)
     */
    void filter_candidates() {
        std::vector<std::string> to_remove;
        
        for (const auto& entry : candidates_) {
            const auto& site = entry.second;
            
            // Filter 1: Minimum evidence
            if (site.left_clip_count + site.right_clip_count < 3 &&
                site.left_discordant + site.right_discordant < 2) {
                to_remove.push_back(entry.first);
                continue;
            }
            
            // Filter 2: AF cutoff (from xTEA's AF conflict check)
            uint32_t total_effective = site.left_clip_count + site.right_clip_count +
                                        site.left_discordant + site.right_discordant;
            if (total_effective > 0) {
                float af = (float)(site.left_clip_consensus + site.right_clip_consensus) / total_effective;
                if (af < min_af_cutoff_) {
                    to_remove.push_back(entry.first);
                    continue;
                }
            }
            
            // Filter 3: Cluster consistency (for two_side only)
            if (site.support_type == "two_side") {
                uint32_t gap = 0;
                if (site.left_cluster_end > 0 && site.right_cluster_start > 0) {
                    gap = std::abs((int)(site.right_cluster_start - site.left_cluster_end));
                }
                
                uint32_t cutoff = clip_cluster_diff_cutoff_;
                if (site.te_type == "SVA" || site.te_type == "Retroposon") {
                    cutoff = sva_clip_cluster_diff_cutoff_;
                }
                
                if (gap > cutoff) {
                    to_remove.push_back(entry.first);
                    continue;
                }
            }
        }
        
        // Remove filtered candidates
        for (const auto& key : to_remove) {
            candidates_.erase(key);
        }
        
        std::cerr << "Filtered " << to_remove.size() << " candidates. " 
                  << candidates_.size() << " remaining." << std::endl;
    }
    
    /**
     * @brief Write output in BED format (compatible with MEGAnE and xTEA)
     * Format: chr, start, end, info, score, strand, extra_fields
     */
    bool write_output_bed(const std::string& output_path) {
        std::ofstream ofs(output_path);
        if (!ofs) {
            std::cerr << "Error: Cannot write to " << output_path << std::endl;
            return false;
        }
        
        // Write header
        ofs << "#chrom\tstart\tend\tinfo\tgenotype\tstrand\tleft_clip\tright_clip\tleft_disc\tright_disc\tpolyA\tsupport_type\n";
        
        for (const auto& entry : candidates_) {
            const auto& site = entry.second;
            
            // BED is 0-based, half-open
            ofs << site.chrom << "\t"
                << site.pos << "\t"
                << (site.pos + 1) << "\t"  // end = start + 1 (point event)
                << site.te_type << "\t"
                << "0/0\t"  // genotype (unknown at this stage)
                << ".\t"   // strand
                << site.left_clip_count << "\t"
                << site.right_clip_count << "\t"
                << site.left_discordant << "\t"
                << site.right_discordant << "\t"
                << (site.left_polyA + site.right_polyA) << "\t"
                << site.support_type << "\n";
        }
        
        ofs.close();
        std::cerr << "Wrote " << candidates_.size() << " candidates to " << output_path << std::endl;
        return true;
    }
    
    /**
     * @brief Get detection statistics
     */
    void get_statistics(uint64_t& total_reads, uint64_t& candidates_found) {
        if (bam_processor_) {
            auto stats = bam_processor_->get_statistics();
            total_reads = stats.total_reads;
        } else {
            total_reads = 0;
        }
        candidates_found = candidates_.size();
    }
};

} // namespace candidate_detector_hpp

#ifndef BUILD_LIB
// Main function (for standalone testing)
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.bam/cram> <output.bed> [reference.fa] [n_threads]" << std::endl;
        return 1;
    }
    
    std::string bam_path = argv[1];
    std::string output_bed = argv[2];
    std::string ref_path = (argc > 3) ? argv[3] : "";
    int n_threads = (argc > 4) ? atoi(argv[4]) : 1;
    
    candidate_detector_hpp::CandidateDetector detector;
    
    if (!detector.initialize(bam_path, ref_path, n_threads)) {
        std::cerr << "Error: Initialization failed" << std::endl;
        return 1;
    }
    
    if (!detector.detect_candidates(output_bed)) {
        std::cerr << "Error: Detection failed" << std::endl;
        return 1;
    }
    
    uint64_t total_reads, candidates_found;
    detector.get_statistics(total_reads, candidates_found);
    
    std::cerr << "\n=== Detection Summary ===" << std::endl;
    std::cerr << "Total reads processed: " << total_reads << std::endl;
    std::cerr << "Candidates found: " << candidates_found << std::endl;
    std::cerr << "This used SINGLE-PASS approach (vs xTEA's chromosome-splitting)" << std::endl;
    
    return 0;
}
#endif // BUILD_LIB
