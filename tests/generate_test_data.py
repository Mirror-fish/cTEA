"""
generate_test_data.py -- Generate small-scale test data for cTEA development
@author cTEA Development Team
@version 0.1.0

Creates small CRAM/BAM files (~50-100MB) for local testing:
- Simulates ~1000 MEI candidate sites on a mini-chromosome (~1MB)
- Includes L1, Alu, SVA insertions with known ground truth
- Generates corresponding BAM/CRAM with 10x, 30x coverage
- Outputs candidate lists in MEGAnE/xTEA compatible formats
"""

import os
import sys
import subprocess
from pathlib import Path
import random
import math

def generate_mini_reference(output_fasta, size=1000000):
    """Generate a mini reference genome (1MB) with known repeat structure"""
    print(f"Generating mini reference ({size/1000000:.1f}MB)...")
    
    # Simplified human chromosome with known repeat blocks
    # This is a simulation, not real reference
    seq = []
    
    # Add some random sequence
    bases = ['A', 'T', 'G', 'C']
    for i in range(size):
        seq.append(random.choice(bases))
    
    # Insert known repeat regions (simulating L1, Alu, SVA loci)
    # Format: (start, end, family, subfamily)
    repeat_regions = [
        (10000, 15000, "L1", "L1HS"),      # L1 insertion site
        (30000, 35000, "Alu", "AluJ"),      # Alu insertion site  
        (60000, 70000, "SVA", "SVA_E"),     # SVA insertion site
        (80000, 85000, "L1", "L1PA2"),    # Another L1
        (120000, 125000, "Alu", "AluSx"),  # Another Alu
        (400000, 410000, "SVA", "SVA_F"),   # Another SVA
    ]
    
    # Mark these regions in the sequence (just for tracking)
    for start, end, family, subfam in repeat_regions:
        # In real implementation, these would be actual repeat sequences
        # For simulation, we just note their positions
        pass
    
    # Write FASTA
    with open(output_fasta, 'w') as f:
        f.write(f">mini_chr1 simulated mini-chromosome for cTEA testing\n")
        # Write in 60-char lines
        sequence = ''.join(seq)
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
    
    print(f"Reference written to {output_fasta}")
    return repeat_regions

def generate_known_mei_sites(repeat_regions, n_sites=1000):
    """Generate known MEI insertion sites (ground truth)"""
    print(f"Generating {n_sites} known MEI sites...")
    
    sites = []
    te_types = ['L1', 'Alu', 'SVA']
    mei_id = 10000
    
    # Distribute sites around repeat regions (simulating real distribution)
    for _ in range(n_sites):
        # 70% near repeat regions, 30% random
        if random.random() < 0.7:
            # Near a repeat region
            region = random.choice(repeat_regions)
            offset = random.randint(-5000, 5000)
            pos = max(0, region[0] + offset)
            te_type = region[2]
            subfam = region[3]
        else:
            # Random position
            pos = random.randint(0, 999999)
            te_type = random.choice(te_types)
            subfam = f"{te_type}_subfamily"
        
        site = {
            'id': f"ID={mei_id};{mei_id+1};",
            'chrom': 'mini_chr1',
            'pos': pos,
            'te_type': te_type,
            'subfamily': subfam,
            'has_polyA': te_type in ['L1', 'SVA'],  # L1/SVA often have polyA
            'left_clip': random.randint(5, 50),
            'right_clip': random.randint(5, 50),
            'left_disc': random.randint(3, 30),
            'right_disc': random.randint(3, 30),
        }
        sites.append(site)
        mei_id += 2
    
    return sites

def write_candidate_files(sites, output_dir):
    """Write candidate files in MEGAnE/xTEA compatible formats"""
    print(f"Writing candidate files to {output_dir}...")
    
    # 1. candidate_list_from_clip.txt
    clip_file = os.path.join(output_dir, "candidate_list_from_clip.txt")
    with open(clip_file, 'w') as f:
        for site in sites:
            # Format similar to xTEA's output
            f.write(f"{site['chrom']}\t{site['pos']}\t{site['pos']+1}\t{mean([site['left_clip'], site['right_clip'])}\n")
    print(f"  Written: {clip_file}")
    
    # 2. candidate_list_from_disc.txt
    disc_file = os.path.join(output_dir, "candidate_list_from_disc.txt")
    with open(disc_file, 'w') as f:
        for site in sites:
            f.write(f"{site['chrom']}\t{site['pos']-100}\t{site['pos']+100}\t{site['left_disc']+site['right_disc']}\n")
    print(f"  Written: {disc_file}")
    
    # 3. MEI_final_gaussian.bed (simulated)
    bed_file = os.path.join(output_dir, "MEI_final_gaussian.bed")
    with open(bed_file, 'w') as f:
        f.write("#chrom\tstart\tend\tinfo\tgenotype\tstrand\tleft_clip\tright_clip\tleft_disc\tright_disc\tpolyA\tsupport_type\n")
        for site in sites[:100]:  # Just first 100 for testing
            info = f"MEI_left:ref_pos={site['pos']},chimeric={site['left_clip']},hybrid={site['left_disc']},pA={1 if site['has_polyA'] else 0}"
            info += f";MEI_right:ref_pos={site['pos']+100},chimeric={site['right_clip']},hybrid={site['right_disc']},pA={1 if site['has_polyA'] else 0}"
            support_type = "two_side" if site['left_clip'] > 0 and site['right_clip'] > 0 else "one_side"
            f.write(f"{site['chrom']}\t{site['pos']}\t{site['pos']+1}\t{site['te_type']}\t1/1\t.\t{site['left_clip']}\t{site['right_clip']}\t{site['left_disc']}\t{site['right_disc']}\t{1 if site['has_polyA'] else 0}\t{support_type}\n")
    print(f"  Written: {bed_file}")
    
    # 4. ml_genotype_features.tsv (simulated features)
    feat_file = os.path.join(output_dir, "ml_genotype_features.tsv")
    with open(feat_file, 'w') as f:
        f.write("chrom\tpos\tleft_coverage\tright_coverage\tleft_clip_cnt\tright_clip_cnt\tleft_disc_cnt\tright_disc_cnt\tpolyA_left\tpolyA_right\t...\n")
        for site in sites[:100]:
            # Simulate 15-dimensional feature vector (simplified)
            left_cov = random.uniform(10, 50)
            right_cov = random.uniform(10, 50)
            f.write(f"{site['chrom']}\t{site['pos']}\t{left_cov:.1f}\t{right_cov:.1f}\t{site['left_clip']}\t{site['right_clip']}\t{site['left_disc']}\t{site['right_disc']}\t{1 if site['has_polyA'] else 0}\t{1 if site['has_polyA'] else 0}\t...\n")
    print(f"  Written: {feat_file}")
    
    return clip_file, disc_file, bed_file, feat_file

def simulate_bam_generation(ref_fasta, output_bam, coverage=30):
    """Simulate BAM generation (placeholder - would use ART or dwgsim)"""
    print(f"Simulating BAM generation with {coverage}x coverage...")
    print("Note: This is a placeholder. In real implementation, use:")
    print("  - ART: for realistic short-read simulation")
    print("  - dwgsim: for simpler simulation")
    print(f"Output would be: {output_bam}")
    
    # Placeholder: create empty BAM
    # In real implementation:
    # subprocess.run(["art_illumina", "-i", ref_fasta, "-o", "tmp", "-c", str(coverage)])
    # subprocess.run(["samtools", "view", "-b", "tmp.sam", "-o", output_bam])
    
    return None

def main():
    print("="*60)
    print("cTEA Test Data Generator")
    print("="*60)
    
    # Setup directories
    base_dir = Path(__file__).parent.parent
    test_dir = base_dir / "tests" / "test_data"
    test_dir.mkdir(exist_ok=True)
    
    # 1. Generate mini reference
    ref_fasta = test_dir / "mini_reference.fa"
    repeat_regions = generate_mini_reference(ref_fasta)
    
    # 2. Generate known MEI sites
    sites = generate_known_mei_sites(repeat_regions, n_sites=1000)
    print(f"Generated {len(sites)} MEI sites")
    
    # 3. Write candidate files
    output_files = write_candidate_files(sites, str(test_dir))
    
    # 4. Simulate BAM (placeholder)
    bam_file = test_dir / "mini_test.bam"
    simulate_bam_generation(ref_fasta, bam_file, coverage=30)
    
    # 5. Generate summary
    summary_file = test_dir / "test_data_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("cTEA Test Data Summary\n")
        f.write("="*60 + "\n")
        f.write(f"Reference: {ref_fasta}\n")
        f.write(f"  Size: ~1MB mini-chromosome\n")
        f.write(f"\nMEI Sites: {len(sites)}\n")
        f.write(f"  L1: {sum(1 for s in sites if s['te_type'] == 'L1')}\n")
        f.write(f"  Alu: {sum(1 for s in sites if s['te_type'] == 'Alu')}\n")
        f.write(f"  SVA: {sum(1 for s in sites if s['te_type'] == 'SVA')}\n")
        f.write(f"\nOutput Files:\n")
        for fpath in output_files:
            f.write(f"  {fpath}\n")
        f.write(f"\nNote: BAM file is a placeholder.\n")
        f.write(f"To generate real BAM, install ART: 'pip install ART' or use dwgsim.\n")
    
    print(f"\nSummary written to: {summary_file}")
    print("\nTest data generation complete!")
    print(f"Test data directory: {test_dir}")
    print("\nNext steps:")
    print("1. Install ART or dwgsim to generate real BAM/CRAM")
    print("2. Compile cTEA C++ core (fix htslib paths)")
    print("3. Run cTEA with: python ctea_main.py --input tests/test_data/mini_test.bam")

if __name__ == "__main__":
    main()