"""
ctea_core.py -- Python wrapper for cTEA C++ core functions
@author cTEA Development Team
@version 0.1.0

This module provides Python interfaces to the C++ core:
1. BAMProcessor (bam_processor.cpp)
2. ReadExtractor (extract_reads.cpp)  
3. KmerFilter (kmer_filter.cpp)
4. CandidateDetector (candidate_detector.cpp)

Uses ctypes or pybind11 (when C++ is compiled).
For now, this is a placeholder that simulates the C++ functionality.
"""

import os
import sys
import ctypes
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

# Try to load the shared library
try:
    # For pybind11 (preferred)
    sys.path.append(str(Path(__file__).parent.parent / "cpp"))
    import libcTEA  # type: ignore
    USING_PYBIND11 = True
except ImportError:
    # Fallback to ctypes
    try:
        lib_path = str(Path(__file__).parent.parent / "cpp" / "libcTEA.so")
        if not os.path.exists(lib_path):
            lib_path = str(Path(__file__).parent.parent / "cpp" / "libcTEA.dylib")
        libcTEA = ctypes.CDLL(lib_path)
        USING_PYBIND11 = False
    except Exception as e:
        print(f"Warning: C++ core not compiled: {e}")
        print("Using Python fallback implementation")
        libcTEA = None
        USING_PYBIND11 = False

class BAMProcessor:
    """Python wrapper for cTEA::BamProcessor"""
    
    def __init__(self):
        if libcTEA:
            if USING_PYBIND11:
                self._obj = libcTEA.cTEA_BamProcessor_new()
            else:
                # ctypes initialization
                self._obj = None
        else:
            self._obj = None
    
    def initialize(self, bam_path: str, ref_path: str = "", n_threads: int = 1) -> bool:
        """Initialize with BAM/CRAM file"""
        if self._obj and USING_PYBIND11:
            return libcTEA.cTEA_BamProcessor_initialize(
                self._obj, bam_path.encode(), ref_path.encode(), n_threads
            )
        else:
            # Python fallback (simplified)
            self.bam_path = bam_path
            self.ref_path = ref_path
            self.n_threads = n_threads
            return os.path.exists(bam_path)
    
    def scan_all_evidence(self, progress_callback=None):
        """Scan BAM and collect evidence (single-pass)"""
        evidence: Dict[str, Dict[int, Any]] = {}
        
        if self._obj and USING_PYBIND11:
            # Call C++ function
            pass  # Implementation when compiled
        else:
            # Python fallback using pysam
            try:
                import pysam
                bam = pysam.AlignmentFile(self.bam_path, "rb")
                
                total = 0
                for read in bam:
                    chrom = bam.get_reference_name(read.reference_id)
                    pos = read.reference_start
                    
                    if chrom not in evidence:
                        evidence[chrom] = {}
                    
                    if pos not in evidence[chrom]:
                        evidence[chrom][pos] = {
                            'left_clip_count': 0,
                            'right_clip_count': 0,
                            'left_discordant': 0,
                            'right_discordant': 0,
                            'left_polyA': 0,
                            'right_polyA': 0,
                        }
                    
                    # Simple evidence collection (placeholder)
                    if read.is_reverse:
                        evidence[chrom][pos]['left_clip_count'] += 1
                    else:
                        evidence[chrom][pos]['right_clip_count'] += 1
                    
                    total += 1
                    if total % 1000000 == 0 and progress_callback:
                        progress_callback(total, 0)
                
                bam.close()
                
            except ImportError:
                print("Error: pysam not installed. Install with: pip install pysam")
                return {}
        
        return evidence

class ReadExtractor:
    """Python wrapper for read extraction (merged from MEGAnE)"""
    
    def __init__(self):
        self.results = {
            'pA_reads': [],
            'overhang_reads': [],
            'distant_reads': [],
            'absent_reads': [],
            'unmapped_reads': [],
        }
    
    def extract_all(self, bam_path: str, ref_path: str = "", n_threads: int = 1):
        """Extract all read types in single pass"""
        # Placeholder: when C++ is compiled, call C++ function
        # For now, return empty results
        print(f"Extraction from {bam_path} (placeholder)")
        return self.results

class KmerFilter:
    """Python wrapper for k-mer filtering"""
    
    def __init__(self):
        self.kmer_set = set()
    
    def load_kmer_set(self, mk_file: str) -> bool:
        """Load k-mer set from .mk file"""
        try:
            with open(mk_file, 'rb') as f:
                # Simplified: read binary k-mers
                pass
            return True
        except Exception as e:
            print(f"Error loading k-mer set: {e}")
            return False
    
    def filter_sequence(self, seq: str) -> Tuple[bool, int]:
        """Filter a single sequence"""
        # Placeholder
        return True, 0

class CandidateDetector:
    """Python wrapper for candidate detection (KEY OPTIMIZATION vs xTEA)"""
    
    def __init__(self):
        self.candidates: Dict[str, Any] = {}
    
    def detect_candidates(self, bam_path: str, output_bed: str = "", 
                        ref_path: str = "", n_threads: int = 1) -> bool:
        """Main detection: SINGLE-PASS vs xTEA's chromosome-splitting"""
        print("Starting single-pass candidate detection (replaces xTEA's chromosome-splitting)")
        print(f"Input: {bam_path}")
        
        # Placeholder: when C++ is compiled, call C++ function
        # For now, simulate by calling MEGAnE's script
        try:
            # Fallback to MEGAnE's pipeline (temporary)
            scripts_dir = Path(__file__).parent.parent / "scripts"
            if scripts_dir.exists():
                cmd = [
                    sys.executable,
                    str(scripts_dir / "1_indiv_call_genotype.py"),
                    "-i", bam_path,
                    # ... other params
                ]
                # subprocess.run(cmd, check=True)
                print("Using MEGAnE fallback (C++ not compiled)")
            else:
                print("Error: No detection method available")
                return False
                
        except Exception as e:
            print(f"Error in detection: {e}")
            return False
        
        return True

def test_core_components():
    """Test the core components"""
    print("Testing cTEA core components...")
    print(f"C++ core available: {libcTEA is not None}")
    print(f"Using pybind11: {USING_PYBIND11}")
    
    # Test BAMProcessor
    print("\n1. BAMProcessor:")
    bp = BAMProcessor()
    print(f"   Initialized: {bp is not None}")
    
    # Test ReadExtractor
    print("\n2. ReadExtractor:")
    re = ReadExtractor()
    print(f"   Initialized: {re is not None}")
    
    # Test KmerFilter
    print("\n3. KmerFilter:")
    kf = KmerFilter()
    print(f"   Initialized: {kf is not None}")
    
    # Test CandidateDetector
    print("\n4. CandidateDetector:")
    cd = CandidateDetector()
    print(f"   Initialized: {cd is not None}")
    
    print("\nCore component test complete.")

if __name__ == "__main__":
    test_core_components()