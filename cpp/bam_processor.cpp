/**
 * @file bam_processor.cpp
 * @brief Implementation of BAM/CRAM processor using htslib
 * @author cTEA Development Team
 * @version 0.1.0
 * 
 * Key optimizations vs xTEA:
 * 1. Uses htslib C API (vs xTEA's Python pysam)
 * 2. Single-pass scanning (vs xTEA's chromosome-split + merge)
 * 3. OpenMP threading (vs xTEA's Python multiprocessing)
 * 4. In-memory hash maps (vs xTEA's temporary files)
 */

#include "bam_processor.h"
#include <htslib/sam.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>

namespace cTEA {

// Mutex for thread-safe operations
static std::mutex g_mutex;

BamProcessor::BamProcessor() 
    : bam_file_(nullptr),
      bam_header_(nullptr),
      bam_index_(nullptr),
      n_threads_(1) {
    stats_ = {0, 0, 0, 0, 0, 0.0, 0.0};
}

BamProcessor::~BamProcessor() {
    close_bam_file();
}

bool BamProcessor::initialize(const std::string& bam_path,
                              const std::string& reference_path,
                              int n_threads) {
    bam_path_ = bam_path;
    reference_path_ = reference_path;
    n_threads_ = n_threads;
    
    return open_bam_file();
}

bool BamProcessor::open_bam_file() {
    // Check file extension to determine if CRAM
    bool is_cram = false;
    size_t dot_pos = bam_path_.find_last_of(".");
    if (dot_pos != std::string::npos) {
        std::string ext = bam_path_.substr(dot_pos);
        if (ext == ".cram" || ext == ".CRAM") {
            is_cram = true;
        }
    }
    
    // Set up htslib thread pool
    // Use OMP_NUM_THREADS environment variable if n_threads_ <= 1
    if (n_threads_ <= 1) {
        const char* env_threads = getenv("OMP_NUM_THREADS");
        if (env_threads != nullptr) {
            n_threads_ = atoi(env_threads);
        }
    }
    
    // Open file
    // For CRAM, we need reference FASTA
    if (is_cram && !reference_path_.empty()) {
        // Set reference for CRAM decoding
        std::string ref_opt = "reference=" + reference_path_;
        // New htslib API: sam_open takes only 2 args, set reference via hts_set_fai_filename
        bam_file_ = sam_open(bam_path_.c_str(), "r");
        if (bam_file_ != nullptr) {
            hts_set_fai_filename(bam_file_, reference_path_.c_str());
        }
    } else {
        bam_file_ = sam_open(bam_path_.c_str(), "r");
    }
    
    if (bam_file_ == nullptr) {
        std::cerr << "Error: Cannot open file " << bam_path_ << std::endl;
        return false;
    }
    
    // Set thread pool if supported
    if (n_threads_ > 1) {
        // Use hts_set_threads instead of hts_set_thread_pool for compatibility
        hts_set_threads(bam_file_, n_threads_);
    }
    
    // Read header
    bam_header_ = sam_hdr_read(bam_file_);
    if (bam_header_ == nullptr) {
        std::cerr << "Error: Cannot read header from " << bam_path_ << std::endl;
        close_bam_file();
        return false;
    }
    
    // Load index if available (for random access)
    bam_index_ = sam_index_load(bam_file_, bam_path_.c_str());
    // Index is optional - we can still do sequential scan
    
    return true;
}

void BamProcessor::close_bam_file() {
    if (bam_index_ != nullptr) {
        hts_idx_destroy(bam_index_);
        bam_index_ = nullptr;
    }
    if (bam_header_ != nullptr) {
        bam_hdr_destroy(bam_header_);
        bam_header_ = nullptr;
    }
    if (bam_file_ != nullptr) {
        sam_close(bam_file_);
        bam_file_ = nullptr;
    }
}

BamHeaderInfo BamProcessor::get_header_info() const {
    BamHeaderInfo info;
    if (bam_header_ == nullptr) {
        return info;
    }
    
    // Parse header
    info.version = "1.0";  // Default version
    info.n_chromosomes = bam_header_->n_targets;
    
    for (int i = 0; i < bam_header_->n_targets; i++) {
        info.chromosome_names.push_back(bam_header_->target_name[i]);
        info.chromosome_lengths.push_back(bam_header_->target_len[i]);
    }
    
    info.reference_fasta = reference_path_;
    return info;
}

bool BamProcessor::scan_all_evidence(
    std::unordered_map<std::string, 
                      std::unordered_map<uint32_t, CandidateEvidence>>& evidence_output,
    std::function<void(uint32_t, uint32_t)> progress_callback) {
    
    if (bam_file_ == nullptr || bam_header_ == nullptr) {
        std::cerr << "Error: BAM file not initialized" << std::endl;
        return false;
    }
    
    // Set the evidence pointer for process_alignment to use
    evidence_ptr_ = &evidence_output;
    
    bam1_t* aln = bam_init1();
    if (aln == nullptr) {
        std::cerr << "Error: Cannot allocate alignment" << std::endl;
        evidence_ptr_ = nullptr;
        return false;
    }
    
    uint64_t read_count = 0;
    uint64_t total_reads = 0;
    int ret;
    
    std::cerr << "Scanning BAM/CRAM and collecting evidence..." << std::endl;
    
    // Single-pass scan - KEY OPTIMIZATION vs xTEA
    // xTEA does: for chrm in chrms: for each BAM: read temp files, filter by chrm
    // cTEA does: single scan, populate hash maps in memory
    
    // Reuse the same alignment object (don't destroy/reinit each time)
    while ((ret = sam_read1(bam_file_, bam_header_, aln)) > 0) {  // Use > 0 instead of >= 0
        read_count++;
        total_reads++;
        
        // Skip unmapped reads (tid < 0)
        if (aln->core.tid < 0) {
            continue;
        }
        
        // Process alignment and collect evidence
        process_alignment(aln, bam_header_);
        
        // Progress reporting (every 1M reads)
        if (read_count >= 1000000) {
            if (progress_callback) {
                progress_callback(total_reads, 0);  // 0 = unknown total
            }
            read_count = 0;
        }
    }
    
    if (ret < -1) {  // -1 is normal EOF, < -1 is error
        std::cerr << "\nError: sam_read1 failed with return code " << ret 
                  << " at read " << total_reads << std::endl;
    }
    
    if (aln != nullptr) {
        bam_destroy1(aln);
    }
    
    // Clear the evidence pointer
    evidence_ptr_ = nullptr;
    
    stats_.total_reads = total_reads;
    std::cerr << "Scan complete. Total reads: " << total_reads << std::endl;
    std::cerr << "Evidence sites collected: " << evidence_output.size() << " chromosomes" << std::endl;
    return true;
}

void BamProcessor::process_alignment(bam1_t* aln, const bam_hdr_t* header) {
    if (aln == nullptr || header == nullptr) return;
    
    // Skip if no evidence pointer set
    if (evidence_ptr_ == nullptr) {
        // Just update statistics without collecting evidence
        stats_.mapped_reads++;
        if (aln->core.flag & BAM_FUNMAP) {
            stats_.unmapped_reads++;
        }
        return;
    }
    
    // Get basic info
    std::string chromosome = header->target_name[aln->core.tid];
    uint32_t pos = aln->core.pos;
    
    // Check for soft clipping
    bool has_clip = false;
    uint32_t clip_len = 0;
    bool is_left_clip = false;
    
    if (aln->core.n_cigar > 0) {
        uint32_t* cigar = bam_get_cigar(aln);
        
        // Check left clip (first CIGAR op)
        uint32_t first_op = bam_cigar_op(cigar[0]);
        if (first_op == BAM_CSOFT_CLIP) {
            clip_len = bam_cigar_oplen(cigar[0]);
            if (clip_len >= 2) {  // Changed from 5 to 2 to match xTEA's MIN_CLIP_FOR_CANDIDATE
                has_clip = true;
                is_left_clip = true;
            }
        }
        
        // Check right clip (last CIGAR op)
        if (!has_clip) {
            uint32_t last_op = bam_cigar_op(cigar[aln->core.n_cigar - 1]);
            if (last_op == BAM_CSOFT_CLIP) {
                clip_len = bam_cigar_oplen(cigar[aln->core.n_cigar - 1]);
                if (clip_len >= 2) {  // Changed from 5 to 2 to match xTEA
                    has_clip = true;
                }
            }
        }
    }
    
    // Check for discordant pair
    bool is_discordant = false;
    if (!(aln->core.flag & BAM_FUNMAP) && !(aln->core.flag & BAM_FMUNMAP)) {
        bool read_reverse = (aln->core.flag & BAM_FREVERSE) != 0;
        bool mate_reverse = (aln->core.flag & BAM_FMREVERSE) != 0;
        
        // Discordant if same orientation
        if (read_reverse == mate_reverse) {
            is_discordant = true;
        }
        
        // Discordant if large insert size
        if (abs(aln->core.isize) > 1000) {
            is_discordant = true;
        }
    } else if (aln->core.flag & BAM_FMUNMAP) {
        // Unmapped mate = discordant
        is_discordant = true;
    }
    
    // Collect evidence if we have any
    if (has_clip || is_discordant) {
        // Use mutex for thread safety (if multi-threaded)
        std::lock_guard<std::mutex> lock(g_mutex);
        
        // Get or create the position map for this chromosome
        auto& pos_map = (*evidence_ptr_)[chromosome];
        
        // Get or create evidence for this position
        auto it = pos_map.find(pos);
        if (it == pos_map.end()) {
            // Create new evidence
            CandidateEvidence ev;
            ev.chromosome = chromosome;
            ev.position = pos;
            ev.left_clip_count = 0;
            ev.right_clip_count = 0;
            ev.left_clip_consensus = 0;
            ev.right_clip_consensus = 0;
            ev.left_discordant = 0;
            ev.right_discordant = 0;
            ev.left_polyA = 0;
            ev.right_polyA = 0;
            ev.left_coverage = 0.0;
            ev.right_coverage = 0.0;
            ev.te_type = "";
            ev.subfamily = "";
            ev.divergence = 100.0;
            pos_map[pos] = ev;
            it = pos_map.find(pos);
        }
        
        // Update evidence
        if (has_clip) {
            if (is_left_clip) {
                it->second.left_clip_count++;
            } else {
                it->second.right_clip_count++;
            }
        }
        
        if (is_discordant) {
            if (is_left_clip || (aln->core.flag & BAM_FREVERSE) == 0) {
                it->second.left_discordant++;
            } else {
                it->second.right_discordant++;
            }
        }
    }
    
    // Update statistics
    stats_.mapped_reads++;
    if (has_clip) {
        stats_.clipped_reads++;
    }
    if (aln->core.flag & BAM_FUNMAP) {
        stats_.unmapped_reads++;
    }
    if (is_discordant) {
        stats_.discordant_pairs++;
    }
}

void BamProcessor::extract_clip_info(bam1_t* aln, AlignmentRecord& record) {
    if (aln->core.n_cigar == 0) return;
    
    uint32_t* cigar = bam_get_cigar(aln);
    int left_clip = 0, right_clip = 0;
    
    // Check left clip (first CIGAR op)
    uint32_t first_op = bam_cigar_op(cigar[0]);
    if (first_op == BAM_CSOFT_CLIP) {
        left_clip = bam_cigar_oplen(cigar[0]);
        record.is_left_clipped = true;
    }
    
    // Check right clip (last CIGAR op)
    uint32_t last_op = bam_cigar_op(cigar[aln->core.n_cigar - 1]);
    if (last_op == BAM_CSOFT_CLIP) {
        right_clip = bam_cigar_oplen(cigar[aln->core.n_cigar - 1]);
        record.is_right_clipped = true;
    }
    
    // Extract clipped sequences (simplified - don't use bam_seq_const)
    if ((left_clip > 0 || right_clip > 0) && aln->core.l_qseq > 0) {
        // Get sequence using bam_get_seq (returns uint8_t*)
        uint8_t* seq_ptr = bam_get_seq(aln);
        if (seq_ptr != nullptr) {
            // For simplicity, just note that clipping exists
            // Full sequence extraction would need seq_nt16_str conversion
            record.clipped_seq_left = left_clip > 0 ? "CLIPPED_LEFT_" + std::to_string(left_clip) : "";
            record.clipped_seq_right = right_clip > 0 ? "CLIPPED_RIGHT_" + std::to_string(right_clip) : "";
        }
    }
}

void BamProcessor::extract_disc_info(bam1_t* aln, AlignmentRecord& record) {
    // Check if discordant pair
    // Discordant = unmapped mate, or wrong orientation, or large insert size
    if (record.mate_is_unmapped) {
        stats_.discordant_pairs++;
        return;
    }
    
    // Check orientation (F2F1 = discordant for paired-end)
    bool read_reverse = (aln->core.flag & BAM_FREVERSE) != 0;
    bool mate_reverse = (aln->core.flag & BAM_FMREVERSE) != 0;
    
    if (read_reverse == mate_reverse) {
        // Same orientation = discordant
        stats_.discordant_pairs++;
        return;
    }
    
    // Check insert size
    if (abs(record.insert_size) > 1000) {  // Threshold for discordant
        stats_.discordant_pairs++;
    }
}

std::string BamProcessor::get_clipped_sequence(bam1_t* aln, bool left_clip) {
    if (aln->core.l_qseq == 0) return "";
    
    uint32_t* cigar = bam_get_cigar(aln);
    int clip_len = 0;
    
    if (left_clip) {
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
            clip_len = bam_cigar_oplen(cigar[0]);
        }
        if (clip_len > 0 && clip_len <= aln->core.l_qseq) {
            // Use bam_get_seq and convert to string
            uint8_t* seq = bam_get_seq(aln);
            return std::string(seq, seq + clip_len);
        }
    } else {
        if (bam_cigar_op(cigar[aln->core.n_cigar - 1]) == BAM_CSOFT_CLIP) {
            clip_len = bam_cigar_oplen(cigar[aln->core.n_cigar - 1]);
        }
        if (clip_len > 0 && clip_len <= aln->core.l_qseq) {
            uint8_t* seq = bam_get_seq(aln);
            return std::string(seq + aln->core.l_qseq - clip_len, seq + aln->core.l_qseq);
        }
    }
    
    return "";
}

bool BamProcessor::extract_clipped_reads(const std::string& output_fasta, int min_clip_length) {
    // Re-implement MEGAnE's extract_discordant.cpp logic
    // but optimized for single-pass + in-memory output
    if (bam_file_ == nullptr) return false;
    
    std::ofstream fasta_out(output_fasta);
    if (!fasta_out.is_open()) {
        std::cerr << "Error: Cannot open " << output_fasta << std::endl;
        return false;
    }
    
    bam1_t* aln = bam_init1();
    int read_count = 0;
    
    while (sam_read1(bam_file_, bam_header_, aln) >= 0) {
        if (aln->core.l_qseq == 0) {
            bam_destroy1(aln);
            aln = bam_init1();
            continue;
        }
        
        // Check for soft-clips
        uint32_t* cigar = bam_get_cigar(aln);
        bool has_clip = false;
        int clip_start = -1, clip_len = 0;
        
        // Check left clip
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
            clip_len = bam_cigar_oplen(cigar[0]);
            if (clip_len >= min_clip_length) {
                has_clip = true;
                clip_start = 0;
            }
        }
        
        // Check right clip
        if (!has_clip && bam_cigar_op(cigar[aln->core.n_cigar - 1]) == BAM_CSOFT_CLIP) {
            clip_len = bam_cigar_oplen(cigar[aln->core.n_cigar - 1]);
            if (clip_len >= min_clip_length) {
                has_clip = true;
                clip_start = aln->core.l_qseq - clip_len;
            }
        }
        
        if (has_clip && clip_len > 0) {
            // Output in FASTA format
            std::string read_name = bam_get_qname(aln);
            uint8_t* seq_ptr = bam_get_seq(aln);
            if (seq_ptr != nullptr && clip_start >= 0 && clip_len > 0) {
                // Convert uint8_t* to std::string (simplified)
                std::string clip_seq;
                for (int i = clip_start; i < clip_start + clip_len && i < aln->core.l_qseq; i++) {
                    clip_seq += "ACGT"[seq_ptr[i] & 0x3];  // Simplified conversion
                }
                
                fasta_out << ">" << read_name << "_" << (clip_start == 0 ? "left" : "right") << "\n";
                // Output in 60-char lines
                for (size_t i = 0; i < clip_seq.length(); i += 60) {
                    fasta_out << clip_seq.substr(i, 60) << "\n";
                }
            }
            read_count++;
        }
        
        bam_destroy1(aln);
        aln = bam_init1();
    }
    
    if (aln != nullptr) bam_destroy1(aln);
    fasta_out.close();
    
    std::cout << "Extracted " << read_count << " clipped reads to " << output_fasta << std::endl;
    return true;
}

// Helper typedef for complex template
typedef std::pair<uint32_t, uint32_t> DiscPair;
typedef std::vector<DiscPair> DiscPairVector;
typedef std::unordered_map<std::string, DiscPairVector> DiscPairMap;

bool BamProcessor::extract_discordant_pairs(
    DiscPairMap& disc_pairs) {
    // Placeholder for discordant pair extraction
    // Will be merged from MEGAnE's extract_discordant.cpp logic
    return true;
}

BamProcessor::ScanStatistics BamProcessor::get_statistics() const {
    return stats_;
}

void BamProcessor::set_region(const std::string& chromosome, uint32_t start, uint32_t end) {
    // For random access if index is available
    if (bam_index_ != nullptr) {
        // htslib supports region string like "chr1:1000-2000"
        // This is a placeholder for future implementation
    }
}

bool BamProcessor::is_cram() const {
    if (bam_file_ == nullptr) return false;
    // Check file format
    // htslib doesn't expose format directly, so check extension
    size_t dot_pos = bam_path_.find_last_of(".");
    if (dot_pos != std::string::npos) {
        std::string ext = bam_path_.substr(dot_pos);
        return (ext == ".cram" || ext == ".CRAM");
    }
    return false;
}

int BamProcessor::get_cram_compression_level() const {
    // CRAM compression level is not directly exposed by htslib
    // This would need to parse the file header or use hts_get_opt()
    return 0;  // Placeholder
}

} // namespace cTEA