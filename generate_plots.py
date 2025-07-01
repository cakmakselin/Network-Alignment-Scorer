#!/usr/bin/env python3
"""
Script to generate plots for the README.
"""

from network_alignment_scorer import NetworkAlignmentScorer, AlignmentPlotter
from pathlib import Path

def main():
    """Generate plots for the alignment results."""
    
    # Initialize scorer and plotter
    scorer = NetworkAlignmentScorer()
    plotter = AlignmentPlotter()
    
    # Score the alignment
    print("Scoring alignment...")
    results = scorer.score_alignment(
        alignment_file=Path("data/alignments/rno-mmu.sif"),
        go_file_1=Path("data/go_annotations/rno.go"),
        go_file_2=Path("data/go_annotations/mmu.go"),
        mapping_file_1=Path("data/mappings/rno.map"),
        mapping_file_2=Path("data/mappings/mmu.map")
    )
    
    # Create output directory for plots
    output_dir = Path("plots")
    output_dir.mkdir(exist_ok=True)
    
    # Generate plots
    print("Generating plots...")
    
    # Coverage breakdown plot
    plotter.plot_coverage_breakdown(results, save_path=output_dir / "coverage_breakdown.png")
    
    # Similarity distribution plot
    plotter.plot_similarity_distribution(results, save_path=output_dir / "similarity_distribution.png")
    
    # Quality metrics plot
    plotter.plot_quality_metrics(results, save_path=output_dir / "quality_metrics.png")
    
    print(f"Plots saved to {output_dir}/")
    print("You can now include these plots in your README.")

if __name__ == "__main__":
    main() 