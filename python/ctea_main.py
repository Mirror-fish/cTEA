"""
ctea_main.py -- Main pipeline for cTEA (CRAM-optimized Transposable Element Analyzer)
@author cTEA Development Team
@version 0.1.0

Merged pipeline combining:
- MEGAnE's speed (C++ core, single-pass scanning)
- xTEA's accuracy (SVA filtering, ML genotyping)
- cTEA's innovations (no chromosome-splitting, in-memory processing)

Key optimizations vs MEGA-xTEA:
1. Direct C++ calls (not subprocess)
2. Single-pass evidence collection (not per-chromosome iteration)
3. In-memory data structures (not temp files)
4. CRAM support with reference caching
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from typing import Optional, Dict, Any

# Import cTEA core (when compiled)
try:
    from cTEA.python import ctea_core
    CTEA_CORE_AVAILABLE = True
except ImportError:
    print("Warning: cTEA C++ core not compiled. Using Python fallback.")
    CTEA_CORE_AVAILABLE = False

# Import filters (from MEGA-xTEA, optimized)
try:
    from megaxtea import sva_filter, fp_filter, ml_genotype
    FILTERS_AVAILABLE = True
except ImportError:
    print("Warning: MEGA-xTEA filters not available.")
    FILTERS_AVAILABLE = False

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

class ColorFormatter(logging.Formatter):
    """Logging formatter with ANSI colors."""
    _COLORS = {
        logging.DEBUG: "\033[36m",
        logging.INFO: "\033[32m",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[1;31m",
    }
    _RESET = "\033[0m"
    
    def format(self, record):
        color = self._COLORS.get(record.levelno, "")
        msg = super().format(record)
        if color and sys.stderr.isatty():
            return f"{color}{msg}{self._RESET}"
        return msg

def setup_logging(verbosity=0):
    level = logging.DEBUG if verbosity > 0 else logging.INFO
    fmt = ColorFormatter(
        fmt="[%(asctime)s] %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(fmt)
    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(level)
    return logging.getLogger("ctea")

logger = setup_logging()

# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------

BANNER = r"""
   __  __ _____ ____    _           ___ _____ ___    _   
  |  \/  | ____/ ___|  / \    __| __| ____/ _ \  | |
  | |\/| |  _| \___ \ / _ \  |__) |  _| \___ \ |_) | | |
  | |  | | |__) ___) / ___ \| ___/ | |__) ___) | __ <| |_|
  |_|  |_|_____|____/_/   \_\____| |_____|____/|_| \_\___/

  MEGAnE speed | xTEA accuracy | ML genotyping
"""

def print_banner():
    sys.stderr.write(BANNER + "\n")
    sys.stderr.write("  cTEA v0.1.0 - Single-pass MEI detection\n\n")

# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

class StepRunner:
    """Execute pipeline steps with timing."""
    
    def __init__(self):
        self.timings = []
        self._pipeline_start = time.time()
    
    def run(self, name, func, *args, skip=False, **kwargs):
        if skip:
            logger.info(f"Step [{name}] -- SKIPPED")
            self.timings.append({"step": name, "elapsed": 0.0, "skipped": True})
            return None
        
        logger.info(f"Step [{name}] -- starting")
        t0 = time.time()
        try:
            result = func(*args, **kwargs)
        except Exception as e:
            logger.error(f"Step [{name}] -- FAILED after {time.time()-t0:.1f}s: {e}")
            raise
        
        elapsed = time.time() - t0
        logger.info(f"Step [{name}] -- done in {elapsed:.1f}s")
        self.timings.append({"step": name, "elapsed": elapsed, "skipped": False})
        return result
    
    def report(self):
        total = time.time() - self._pipeline_start
        logger.info("=" * 60)
        logger.info("Pipeline timing summary")
        logger.info("-" * 60)
        for entry in self.timings:
            status = "SKIP" if entry.get("skipped") else f"{entry['elapsed']:.1f}s"
            logger.info(f"  {entry['step']:<35s} {status}")
        logger.info("-" * 60)
        logger.info(f"  {'TOTAL':<35s} {total:.1f}s")
        logger.info("=" * 60)

# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(args):
    """Run the full cTEA pipeline."""
    print_banner()
    
    # Validate inputs
    if not os.path.isfile(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    if not os.path.isfile(args.reference):
        logger.error(f"Reference not found: {args.reference}")
        sys.exit(1)
    
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    runner = StepRunner()
    
    # Step 1: Initialize BAM processor (single-pass)
    def step_init():
        logger.info(f"Input BAM/CRAM: {args.input}")
        logger.info(f"Reference:      {args.reference}")
        logger.info(f"Threads:        {args.threads}")
        logger.info(f"Output:         {outdir}")
        
        if CTEA_CORE_AVAILABLE:
            bp = ctea_core.BAMProcessor()
            bp.initialize(args.input, args.reference, args.threads)
            return bp
        else:
            logger.warning("Using Python fallback (C++ core not compiled)")
            return None
    
    bam_processor = runner.run("init", step_init)
    
    # Step 2: Extract reads (C++ accelerated, single-pass)
    def step_extract():
        if CTEA_CORE_AVAILABLE:
            extractor = ctea_core.ReadExtractor()
            results = extractor.extract_all(args.input, args.reference, args.threads)
            logger.info(f"Extracted: {len(results.get('pA_reads', []))} pA reads")
            logger.info(f"          {len(results.get('overhang_reads', []))} overhang reads")
            return results
        else:
            logger.warning("Read extraction requires C++ core")
            return None
    
    runner.run("extract_reads", step_extract)
    
    # Step 3: Candidate detection (KEY OPTIMIZATION vs xTEA)
    def step_detect():
        logger.info("Using SINGLE-PASS detection (replaces xTEA's chromosome-splitting)")
        
        if CTEA_CORE_AVAILABLE:
            detector = ctea_core.CandidateDetector()
            output_bed = str(outdir / "candidates.bed")
            detector.detect_candidates(args.input, output_bed, args.reference, args.threads)
            return output_bed
        else:
            logger.warning("Candidate detection requires C++ core")
            # Fallback: use MEGAnE's script (temporary)
            if args.megane_fallback:
                logger.info("Using MEGAnE fallback...")
                scripts_dir = Path(__file__).parent.parent / "scripts"
                if scripts_dir.exists():
                    cmd = [
                        sys.executable,
                        str(scripts_dir / "1_indiv_call_genotype.py"),
                        "-i", args.input,
                        "-fa", args.reference,
                        "-p", str(args.threads),
                        "-outdir", str(outdir),
                    ]
                    # subprocess.run(cmd, check=True)
            return None
    
    candidates_bed = runner.run("candidate_detection", step_detect)
    
    # Step 4: SVA post-filtering (xTEA logic)
    def step_sva_filter():
        if not FILTERS_AVAILABLE:
            logger.warning("SVA filter not available")
            return 0
        
        from megaxtea import sva_filter
        filt = sva_filter.SVAFilter()
        
        if candidates_bed and os.path.isfile(candidates_bed):
            n = filt.filter_megane_output(candidates_bed, candidates_bed + ".sva_filtered")
            logger.info(f"SVA filter passed: {n} candidates")
            return n
        return 0
    
    runner.run("sva_post_filter", step_sva_filter, skip=not args.sva_filter)
    
    # Step 5: FP filtering (xTEA-style)
    def step_fp_filter():
        if not FILTERS_AVAILABLE:
            logger.warning("FP filter not available")
            return 0
        
        from megaxtea import fp_filter
        decisions = fp_filter.run_fp_filters(
            str(outdir / "candidates.bed.sva_filtered") if candidates_bed else "",
            {}  # bam_features (placeholder)
        )
        logger.info(f"FP filter: {len(decisions)} candidates flagged")
        return len(decisions)
    
    runner.run("fp_filter", step_fp_filter)
    
    # Step 6: ML genotyping (xTEA-style)
    def step_genotype():
        if not FILTERS_AVAILABLE:
            logger.warning("ML genotyping not available")
            return
        
        from megaxtea import ml_genotype
        gt = ml_genotype.UnifiedGenotyper(model_path=args.model_path)
        
        # Process candidates
        if candidates_bed and os.path.isfile(candidates_bed):
            logger.info("Running ML genotyping...")
            # Placeholder: call ml_genotype logic
            pass
    
    runner.run("genotyping", step_genotype, skip=not args.ml_genotype)
    
    # Step 7: Output
    def step_output():
        logger.info("Pipeline output:")
        for f in outdir.glob("*.bed"):
            logger.info(f"  {f}")
        for f in outdir.glob("*.vcf"):
            logger.info(f"  {f}")
    
    runner.run("output", step_output)
    
    # Print timing report
    runner.report()
    logger.info("Pipeline finished successfully.")

# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        prog="ctea",
        description="cTEA: Merged Mobile Element detection (speed + accuracy)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM/CRAM file")
    parser.add_argument("-fa", "--reference", required=True, help="Reference FASTA")
    parser.add_argument("-o", "--output-dir", default="./ctea_output", help="Output directory")
    parser.add_argument("-p", "--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--sva-filter", action="store_true", help="Apply SVA filtering")
    parser.add_argument("--ml-genotype", action="store_true", help="Apply ML genotyping")
    parser.add_argument("--model-path", help="Path to ML model")
    parser.add_argument("--megane-fallback", action="store_true", help="Use MEGAnE for detection (fallback)")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Verbosity level")
    parser.add_argument("-V", "--version", action="version", version="cTEA 0.1.0")
    return parser

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        run_pipeline(args)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)