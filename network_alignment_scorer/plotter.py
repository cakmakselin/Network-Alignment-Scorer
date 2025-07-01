"""
Simple plotting module for network alignment results.

This script provides visualization tools for interpreting the results of network alignment scoring based on Gene Ontology (GO) annotations. It was originally developed as an extension to a university bioinformatics course assignment on network alignment validation.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List
import pandas as pd


class AlignmentPlotter:
    """Simple plotter for alignment results."""
    
    def __init__(self):
        """Initialize the plotter."""
        plt.style.use('default')
        self.colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
    
    def plot_coverage_breakdown(self, results: Dict, save_path: str = None):
        """
        Plot coverage breakdown showing mapped vs unmapped proteins.
        
        Args:
            results: Results dictionary from scorer
            save_path: Path to save the plot
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Coverage pie chart
        labels = ['Successfully Scored', 'Unmappable Species 1', 'Unmappable Species 2']
        sizes = [results['scored_pairs'], results['unmappable_1'], results['unmappable_2']]
        
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=self.colors)
        ax1.set_title('Alignment Coverage Breakdown')
        
        # Bar chart of counts
        x = ['Scored', 'Unmappable 1', 'Unmappable 2']
        y = [results['scored_pairs'], results['unmappable_1'], results['unmappable_2']]
        
        bars = ax2.bar(x, y, color=self.colors)
        ax2.set_title('Protein Pair Counts')
        ax2.set_ylabel('Number of Pairs')
        
        # Add value labels on bars
        for bar, value in zip(bars, y):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01*max(y), 
                    f'{value:,}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show()
    
    def plot_similarity_distribution(self, results: Dict, save_path: str = None):
        """
        Plot similarity score distribution.
        
        Args:
            results: Results dictionary from scorer
            save_path: Path to save the plot
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Mean similarity bar
        metrics = ['Mean Similarity']
        values = [results['mean_score']]
        
        bars = ax1.bar(metrics, values, color=self.colors[0])
        ax1.set_title('Mean Jaccard Similarity')
        ax1.set_ylabel('Similarity Score')
        ax1.set_ylim(0, 1)
        
        # Add value label
        for bar, value in zip(bars, values):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                    f'{value:.4f}', ha='center', va='bottom')
        
        # Coverage bar
        coverage_metrics = ['Coverage']
        coverage_values = [results['coverage']]
        
        bars2 = ax2.bar(coverage_metrics, coverage_values, color=self.colors[1])
        ax2.set_title('Alignment Coverage')
        ax2.set_ylabel('Coverage Percentage')
        ax2.set_ylim(0, 1)
        
        # Add value label
        for bar, value in zip(bars2, coverage_values):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                    f'{value:.2%}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show()
    
    def plot_quality_metrics(self, results: Dict, save_path: str = None):
        """
        Plot comprehensive quality metrics.
        
        Args:
            results: Results dictionary from scorer
            save_path: Path to save the plot
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Total pairs breakdown
        labels = ['Scored', 'Unmappable 1', 'Unmappable 2']
        sizes = [results['scored_pairs'], results['unmappable_1'], results['unmappable_2']]
        
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=self.colors)
        ax1.set_title('Total Alignment Pairs Breakdown')
        
        # Similarity score
        ax2.bar(['Mean Similarity'], [results['mean_score']], color=self.colors[0])
        ax2.set_ylabel('Jaccard Similarity')
        ax2.set_ylim(0, 1)
        ax2.set_title('Functional Similarity')
        ax2.text(0, results['mean_score'] + 0.01, f'{results["mean_score"]:.4f}', 
                ha='center', va='bottom')
        
        # Coverage
        ax3.bar(['Coverage'], [results['coverage']], color=self.colors[1])
        ax3.set_ylabel('Coverage')
        ax3.set_ylim(0, 1)
        ax3.set_title('Alignment Coverage')
        ax3.text(0, results['coverage'] + 0.01, f'{results["coverage"]:.2%}', 
                ha='center', va='bottom')
        
        # Total score
        ax4.bar(['Total Score'], [results['total_score']], color=self.colors[2])
        ax4.set_ylabel('Score')
        ax4.set_title('Cumulative Similarity Score')
        ax4.text(0, results['total_score'] + 0.01*results['total_score'], 
                f'{results["total_score"]:.2f}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show() 