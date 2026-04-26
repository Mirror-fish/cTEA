/**
 * @file bindings.cpp
 * @brief pybind11 bindings for cTEA C++ core
 * @author cTEA Development Team
 * @version 0.1.0
 * 
 * This file exposes cTEA C++ classes to Python via pybind11.
 * Key classes: BamProcessor, CandidateDetector
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <memory>
#include "bam_processor.h"
#include "candidate_detector.h"

namespace py = pybind11;
using namespace cTEA;

PYBIND11_MODULE(libcTEA, m) {
    m.doc() = "cTEA C++ core module - CRAM-optimized Transposable Element Analyzer";
    m.attr("__version__") = "0.1.0";

    // -------------------------------------------------------------------------
    // BamHeaderInfo struct
    // -------------------------------------------------------------------------
    py::class_<BamHeaderInfo>(m, "BamHeaderInfo")
        .def(py::init<>())
        .def_readwrite("version", &BamHeaderInfo::version)
        .def_readwrite("chromosome_names", &BamHeaderInfo::chromosome_names)
        .def_readwrite("chromosome_lengths", &BamHeaderInfo::chromosome_lengths)
        .def_readwrite("n_chromosomes", &BamHeaderInfo::n_chromosomes)
        .def_readwrite("reference_fasta", &BamHeaderInfo::reference_fasta)
        .def("__repr__", [](const BamHeaderInfo &self) {
            return "<BamHeaderInfo: " + std::to_string(self.n_chromosomes) + " chromosomes>";
        });

    // -------------------------------------------------------------------------
    // AlignmentRecord struct
    // -------------------------------------------------------------------------
    py::class_<AlignmentRecord>(m, "AlignmentRecord")
        .def(py::init<>())
        .def_readwrite("chromosome", &AlignmentRecord::chromosome)
        .def_readwrite("pos", &AlignmentRecord::pos)
        .def_readwrite("end", &AlignmentRecord::end)
        .def_readwrite("flag", &AlignmentRecord::flag)
        .def_readwrite("tid", &AlignmentRecord::tid)
        .def_readwrite("mapq", &AlignmentRecord::mapq)
        .def_readwrite("cigar", &AlignmentRecord::cigar)
        .def_readwrite("sequence", &AlignmentRecord::sequence)
        .def_readwrite("qualities", &AlignmentRecord::qualities)
        .def_readwrite("clipped_seq_left", &AlignmentRecord::clipped_seq_left)
        .def_readwrite("clipped_seq_right", &AlignmentRecord::clipped_seq_right)
        .def_readwrite("clip_start", &AlignmentRecord::clip_start)
        .def_readwrite("clip_end", &AlignmentRecord::clip_end)
        .def_readwrite("is_left_clipped", &AlignmentRecord::is_left_clipped)
        .def_readwrite("is_right_clipped", &AlignmentRecord::is_right_clipped)
        .def_readwrite("mate_tid", &AlignmentRecord::mate_tid)
        .def_readwrite("mate_pos", &AlignmentRecord::mate_pos)
        .def_readwrite("insert_size", &AlignmentRecord::insert_size)
        .def_readwrite("is_unmapped", &AlignmentRecord::is_unmapped)
        .def_readwrite("mate_is_unmapped", &AlignmentRecord::mate_is_unmapped)
        .def_readwrite("is_reverse", &AlignmentRecord::is_reverse)
        .def_readwrite("mate_is_reverse", &AlignmentRecord::mate_is_reverse)
        .def_readwrite("tags", &AlignmentRecord::tags);

    // -------------------------------------------------------------------------
    // CandidateEvidence struct
    // -------------------------------------------------------------------------
    py::class_<CandidateEvidence>(m, "CandidateEvidence")
        .def(py::init<>())
        .def_readwrite("chromosome", &CandidateEvidence::chromosome)
        .def_readwrite("position", &CandidateEvidence::position)
        .def_readwrite("left_clip_count", &CandidateEvidence::left_clip_count)
        .def_readwrite("right_clip_count", &CandidateEvidence::right_clip_count)
        .def_readwrite("left_clip_consensus", &CandidateEvidence::left_clip_consensus)
        .def_readwrite("right_clip_consensus", &CandidateEvidence::right_clip_consensus)
        .def_readwrite("left_discordant", &CandidateEvidence::left_discordant)
        .def_readwrite("right_discordant", &CandidateEvidence::right_discordant)
        .def_readwrite("left_polyA", &CandidateEvidence::left_polyA)
        .def_readwrite("right_polyA", &CandidateEvidence::right_polyA)
        .def_readwrite("left_coverage", &CandidateEvidence::left_coverage)
        .def_readwrite("right_coverage", &CandidateEvidence::right_coverage)
        .def_readwrite("te_type", &CandidateEvidence::te_type)
        .def_readwrite("subfamily", &CandidateEvidence::subfamily)
        .def_readwrite("divergence", &CandidateEvidence::divergence)
        .def("__repr__", [](const CandidateEvidence &self) {
            return "<CandidateEvidence: " + self.chromosome + ":" + 
                   std::to_string(self.position) + " " + self.te_type + ">";
        });

    // -------------------------------------------------------------------------
    // BamProcessor::ScanStatistics struct
    // -------------------------------------------------------------------------
    py::class_<BamProcessor::ScanStatistics>(m, "ScanStatistics")
        .def(py::init<>())
        .def_readwrite("total_reads", &BamProcessor::ScanStatistics::total_reads)
        .def_readwrite("mapped_reads", &BamProcessor::ScanStatistics::mapped_reads)
        .def_readwrite("clipped_reads", &BamProcessor::ScanStatistics::clipped_reads)
        .def_readwrite("discordant_pairs", &BamProcessor::ScanStatistics::discordant_pairs)
        .def_readwrite("unmapped_reads", &BamProcessor::ScanStatistics::unmapped_reads)
        .def_readwrite("avg_read_length", &BamProcessor::ScanStatistics::avg_read_length)
        .def_readwrite("estimated_coverage", &BamProcessor::ScanStatistics::estimated_coverage);

    // -------------------------------------------------------------------------
    // BamProcessor class
    // -------------------------------------------------------------------------
    py::class_<BamProcessor>(m, "BamProcessor")
        .def(py::init<>())
        .def("initialize", &BamProcessor::initialize,
             py::arg("bam_path"),
             py::arg("reference_path") = "",
             py::arg("n_threads") = 1,
             R"doc(
             Initialize processor with input file.
             
             Args:
                 bam_path: Path to BAM/CRAM file
                 reference_path: Path to reference FASTA (required for CRAM)
                 n_threads: Number of threads for reading
             
             Returns:
                 True if initialization successful
             )doc")
        .def("get_header_info", &BamProcessor::get_header_info,
             "Get BAM/CRAM header information")
        .def("scan_all_evidence", 
             [](BamProcessor &self, py::function callback) -> py::dict {
                 // Wrap Python callback into C++ std::function
                 std::function<void(uint32_t, uint32_t)> cpp_callback = nullptr;
                 if (callback && !callback.is_none()) {
                     cpp_callback = [&callback](uint32_t current, uint32_t total) {
                         py::gil_scoped_acquire gil;
                         callback(current, total);
                     };
                 }
                 
                 std::unordered_map<std::string, 
                     std::unordered_map<uint32_t, CandidateEvidence>> evidence;
                 
                 // Remove unused variable warning by casting to void
                 (void)self.scan_all_evidence(evidence, cpp_callback);
                 
                 // Convert to Python dict
                 py::dict result;
                 for (const auto &[chrom, pos_map] : evidence) {
                     py::dict chrom_dict;
                     for (const auto &[pos, ev] : pos_map) {
                         py::dict ev_dict;
                         ev_dict["chromosome"] = ev.chromosome;
                         ev_dict["position"] = ev.position;
                         ev_dict["left_clip_count"] = ev.left_clip_count;
                         ev_dict["right_clip_count"] = ev.right_clip_count;
                         ev_dict["left_discordant"] = ev.left_discordant;
                         ev_dict["right_discordant"] = ev.right_discordant;
                         ev_dict["left_polyA"] = ev.left_polyA;
                         ev_dict["right_polyA"] = ev.right_polyA;
                         ev_dict["te_type"] = ev.te_type;
                         ev_dict["subfamily"] = ev.subfamily;
                         ev_dict["divergence"] = ev.divergence;
                         chrom_dict[py::cast(pos)] = ev_dict;
                     }
                     result[py::cast(chrom)] = chrom_dict;
                 }
                 
                 return result;
             },
             py::arg("callback") = py::none(),
             R"doc(
             Single-pass scan to collect all evidence (KEY OPTIMIZATION vs xTEA).
             
             Args:
                 callback: Optional progress callback function(current, total)
             
             Returns:
                 Dictionary of evidence: {chrom: {pos: {evidence_dict}}}
             )doc")
        .def("get_statistics", &BamProcessor::get_statistics,
             "Get statistics from the scan")
        .def("is_cram", &BamProcessor::is_cram,
             "Check if file is CRAM")
        .def("extract_clipped_reads", &BamProcessor::extract_clipped_reads,
             py::arg("output_fasta"),
             py::arg("min_clip_length") = 10)
        .def("extract_discordant_pairs", 
             [](BamProcessor &self) -> py::dict {
                 std::unordered_map<std::string, 
                     std::vector<std::pair<uint32_t, uint32_t>>> disc_pairs;
                 self.extract_discordant_pairs(disc_pairs);
                 
                 py::dict result;
                 for (const auto &[chrom, pairs] : disc_pairs) {
                     py::list pair_list;
                     for (const auto &[start, end] : pairs) {
                         pair_list.append(py::make_tuple(start, end));
                     }
                     result[py::cast(chrom)] = pair_list;
                 }
                 return result;
             },
             "Extract discordant pairs");

    // -------------------------------------------------------------------------
    // CandidateDetector class
    // -------------------------------------------------------------------------
    py::class_<CandidateDetector>(m, "CandidateDetector")
        .def(py::init<>())
        .def("initialize", &CandidateDetector::initialize,
             py::arg("bam_path"),
             py::arg("ref_path") = "",
             py::arg("n_threads") = 1,
             "Initialize detector with BAM/CRAM file")
        .def("detect_candidates", 
             [](CandidateDetector &self, const std::string &output_bed,
                const std::string &rep_kmer_file,
                const std::string &repeat_masker_bed) -> bool {
                 return self.detect_candidates(output_bed, rep_kmer_file, repeat_masker_bed);
             },
             py::arg("output_bed"),
             py::arg("rep_kmer_file") = "",
             py::arg("repeat_masker_bed") = "",
             "Detect TE insertion candidates (single-pass)")
        .def("get_statistics", 
             [](CandidateDetector &self) -> py::dict {
                 uint64_t total_reads, candidates_found;
                 self.get_statistics(total_reads, candidates_found);
                 py::dict result;
                 result["total_reads"] = total_reads;
                 result["candidates_found"] = candidates_found;
                 return result;
             },
             "Get detection statistics");

    // -------------------------------------------------------------------------
    // Utility functions
    // -------------------------------------------------------------------------
    m.def("version", []() -> std::string {
        return "cTEA 0.1.0 - CRAM-optimized Transposable Element Analyzer";
    }, "Get cTEA version string");
}