/**
 * @file bam_processor.h
 * @brief BAM/CRAM processor using htslib
 * @author cTEA Development Team
 * @version 0.1.0
 * 
 * Optimized BAM/CRAM reader with:
 * - htslib backend (supports BAM, CRAM, SAM)
 * - Multi-thread reading via OMP_NUM_THREADS
 * - Memory caching for frequently accessed regions
 * - Avoids xTEA's chromosome-splitting bottleneck
 */

#ifndef CTEA_BAM_PROCESSOR_H
#define CTEA_BAM_PROCESSOR_H

#include <htslib/sam.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <functional>

namespace cTEA {

/**
 * @brief Structure to hold BAM/CRAM header information
 */
struct BamHeaderInfo {
    std::string version;
    std::vector<std::string> chromosome_names;
    std::vector<uint32_t> chromosome_lengths;
    uint32_t n_chromosomes;
    
    // Reference info (needed for CRAM)
    std::string reference_fasta;
};

/**
 * @brief Structure for a single BAM alignment record
 * Optimized for memory efficiency
 */
struct AlignmentRecord {
    std::string chromosome;
    uint32_t pos;              // 0-based position
    uint32_t end;              // calculated end position
    uint16_t flag;
    int32_t tid;                // target id
    uint32_t mapq;
    std::string cigar;
    std::string sequence;
    std::string qualities;
    
    // Clipped sequence info
    std::string clipped_seq_left;
    std::string clipped_seq_right;
    uint32_t clip_start;
    uint32_t clip_end;
    bool is_left_clipped;
    bool is_right_clipped;
    
    // Mate info
    int32_t mate_tid;
    uint32_t mate_pos;
    int32_t insert_size;
    
    // Reference dependent info
    bool is_unmapped;
    bool mate_is_unmapped;
    bool is_reverse;
    bool mate_is_reverse;
    
    // Auxiliary tags
    std::unordered_map<std::string, std::string> tags;
};

/**
 * @brief Candidate site evidence structure
 * Avoids xTEA's multiple temporary files by keeping in memory
 */
struct CandidateEvidence {
    std::string chromosome;
    uint32_t position;
    
    // Clip evidence
    uint32_t left_clip_count;
    uint32_t right_clip_count;
    uint32_t left_clip_consensus;  // clipped part aligned to repeat consensus
    uint32_t right_clip_consensus;
    
    // Discordant pair evidence
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
};

/**
 * @brief Main BAM/CRAM processor class
 * 
 * Key optimizations vs xTEA:
 * 1. Single-pass scanning (vs xTEA's multiple passes)
 * 2. In-memory data structures (vs xTEA's chromosome-splitting + merge)
 * 3. htslib multi-thread (vs xTEA's Python pysam)
 * 4. C++ performance (vs xTEA's Python loops)
 */
class BamProcessor {
public:
    BamProcessor();
    ~BamProcessor();
    
    /**
     * @brief Initialize processor with input file
     * @param bam_path Path to BAM/CRAM file
     * @param reference_path Path to reference FASTA (required for CRAM)
     * @param n_threads Number of threads (from OMP_NUM_THREADS or command line)
     * @return true if initialization successful
     */
    bool initialize(const std::string& bam_path, 
                   const std::string& reference_path = "",
                   int n_threads = 1);
    
    /**
     * @brief Get header information
     */
    BamHeaderInfo get_header_info() const;
    
    /**
     * @brief Single-pass scan to collect all evidence
     * 
     * This is the KEY optimization vs xTEA:
     * xTEA does: scan→chromosome split→process each→merge (x_TEI_locator.py:112-168)
     * cTEA does: single scan→populate hash maps→output (no intermediate files)
     * 
     * @param evidence_output Output map: chromosome → (position → CandidateEvidence)
     * @param callback Optional callback for progress reporting
     */
    bool scan_all_evidence(
        std::unordered_map<std::string, 
                       std::unordered_map<uint32_t, CandidateEvidence>>& evidence_output,
        std::function<void(uint32_t, uint32_t)> progress_callback = nullptr);
    
    /**
     * @brief Extract clipped reads (optimized from MEGAnE's extract_discordant.cpp)
     */
    bool extract_clipped_reads(
        const std::string& output_fasta,
        int min_clip_length = 10);
    
    /**
     * @brief Extract discordant pairs (merged from MEGAnE's extract_discordant + extract_unmapped)
     */
    bool extract_discordant_pairs(
        std::unordered_map<std::string, 
                       std::vector<std::pair<uint32_t, uint32_t>>>& disc_pairs);
    
    /**
     * @brief Get statistics from the scan
     */
    struct ScanStatistics {
        uint64_t total_reads;
        uint64_t mapped_reads;
        uint64_t clipped_reads;
        uint64_t discordant_pairs;
        uint64_t unmapped_reads;
        double avg_read_length;
        double estimated_coverage;
    };
    ScanStatistics get_statistics() const;
    
    /**
     * @brief Set region of interest (for random access if indexed)
     */
    void set_region(const std::string& chromosome, uint32_t start, uint32_t end);
    
    /**
     * @brief Check if file is CRAM
     */
    bool is_cram() const;
    
    /**
     * @brief Get CRAM compression level (for performance tuning)
     */
    int get_cram_compression_level() const;

private:
    // htslib objects
    samFile* bam_file_;
    bam_hdr_t* bam_header_;
    hts_idx_t* bam_index_;
    
    // Configuration
    std::string bam_path_;
    std::string reference_path_;
    int n_threads_;
    
    // Statistics
    ScanStatistics stats_;
    
    // Reference sequence cache (avoid repeated loading)
    std::unordered_map<std::string, std::string> ref_cache_;
    
    // Evidence collection pointer (set by scan_all_evidence)
    std::unordered_map<std::string, 
                      std::unordered_map<uint32_t, CandidateEvidence>>* evidence_ptr_;
    
    // Private methods
    bool open_bam_file();
    void close_bam_file();
    void process_alignment(bam1_t* aln, const bam_hdr_t* header);
    void extract_clip_info(bam1_t* aln, AlignmentRecord& record);
    void extract_disc_info(bam1_t* aln, AlignmentRecord& record);
    std::string get_clipped_sequence(bam1_t* aln, bool left_clip);
    
    // Disabled copy constructor and assignment
    BamProcessor(const BamProcessor&) = delete;
    BamProcessor& operator=(const BamProcessor&) = delete;
};

} // namespace cTEA

#endif // CTEA_BAM_PROCESSOR_H