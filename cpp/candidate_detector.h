/**
 * @file candidate_detector.h
 * @brief Header for CandidateDetector class
 * @author cTEA Development Team
 * @version 0.1.0
 */

#ifndef CTEA_CANDIDATE_DETECTOR_H
#define CTEA_CANDIDATE_DETECTOR_H

#include <string>
#include <unordered_map>
#include <cstdint>
#include "bam_processor.h"

namespace cTEA {

/**
 * @brief Candidate site with accumulated evidence
 */
struct CandidateSite {
    std::string chrom;
    uint32_t pos;
    
    // Evidence counts
    uint32_t left_clip_count;
    uint32_t right_clip_count;
    uint32_t left_clip_consensus;
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
    std::string te_type;
    std::string subfamily;
    float divergence;
    
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
    
    // Support type
    std::string support_type;
    
    CandidateSite();
};

/**
 * @brief Candidate detector class - KEY OPTIMIZATION vs xTEA
 * 
 * Unlike xTEA which splits by chromosome and scans BAM multiple times,
 * this class performs single-pass detection for 3-5x speedup.
 */
class CandidateDetector {
private:
    // In-memory storage
    std::unordered_map<std::string, CandidateSite> candidates_;
    
    // BAM processor
    cTEA::BamProcessor* bam_processor_;
    
    // Parameters
    uint32_t min_clip_len_;
    uint32_t min_disc_count_;
    float min_af_cutoff_;
    uint32_t clip_cluster_diff_cutoff_;
    uint32_t sva_clip_cluster_diff_cutoff_;
    
    // Statistics
    uint64_t total_reads_processed_;
    uint64_t candidates_found_;
    
public:
    CandidateDetector();
    ~CandidateDetector();
    
    bool initialize(const std::string& bam_path,
                  const std::string& ref_path,
                  int n_threads);
    
    bool detect_candidates(const std::string& output_bed,
                         const std::string& rep_kmer_file = "",
                         const std::string& repeat_masker_bed = "");
    
    void get_statistics(uint64_t& total_reads, uint64_t& candidates_found);
    
private:
    void build_candidates(
        const std::unordered_map<std::string, 
                  std::unordered_map<uint32_t, cTEA::CandidateEvidence>>& evidence);
    
    void determine_support_type(CandidateSite& site);
    
    // NEW: Cluster nearby candidates (Priority 1 improvement)
    // Reference: xTEA's chain_regions function (PEAK_WINDOW=100bp)
    void cluster_nearby_candidates();

    // Helper function to process a cluster of nearby candidates
    void process_cluster(
        const std::vector<std::string>& cluster,
        std::unordered_map<std::string, CandidateSite>& clustered_candidates);
    
    void filter_candidates();
    
    bool write_output_bed(const std::string& output_path);
};

} // namespace cTEA

#endif // CTEA_CANDIDATE_DETECTOR_H
