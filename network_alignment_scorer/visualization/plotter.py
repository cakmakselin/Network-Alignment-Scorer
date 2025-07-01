"""
Visualization Module

Provides plotting and visualization capabilities for network alignment results.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

logger = logging.getLogger(__name__)


class AlignmentVisualizer:
    """
    Visualization class for network alignment results.
    
    Provides various plotting capabilities for analyzing and presenting
    alignment quality and statistics.
    """
    
    def __init__(self, style: str = "seaborn-v0_8"):
        """
        Initialize the visualizer.
        
        Args:
            style: Matplotlib style to use
        """
        self.style = style
        plt.style.use(style)
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72',
            'accent': '#F18F01',
            'success': '#C73E1D',
            'neutral': '#6C757D'
        }
        
    def plot_score_distribution(self, scores: List[float], 
                              title: str = "Alignment Score Distribution",
                              save_path: Optional[Path] = None) -> None:
        """
        Plot distribution of alignment scores.
        
        Args:
            scores: List of similarity scores
            title: Plot title
            save_path: Path to save the plot
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Histogram
        ax1.hist(scores, bins=30, alpha=0.7, color=self.colors['primary'], edgecolor='black')
        ax1.set_xlabel('Similarity Score')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Score Distribution Histogram')
        ax1.grid(True, alpha=0.3)
        
        # Box plot
        ax2.boxplot(scores, patch_artist=True, 
                   boxprops=dict(facecolor=self.colors['primary'], alpha=0.7))
        ax2.set_ylabel('Similarity Score')
        ax2.set_title('Score Distribution Box Plot')
        ax2.grid(True, alpha=0.3)
        
        plt.suptitle(title, fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Score distribution plot saved to {save_path}")
        
        plt.show()
    
    def plot_species_comparison(self, comparison_data: Dict[str, any],
                               save_path: Optional[Path] = None) -> None:
        """
        Plot species comparison statistics.
        
        Args:
            comparison_data: Species comparison data
            save_path: Path to save the plot
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Protein counts
        species = ['Species 1', 'Species 2']
        protein_counts = [comparison_data['species1_total_proteins'], 
                         comparison_data['species2_total_proteins']]
        
        bars1 = ax1.bar(species, protein_counts, color=[self.colors['primary'], self.colors['secondary']])
        ax1.set_ylabel('Number of Proteins')
        ax1.set_title('Protein Counts by Species')
        ax1.bar_label(bars1, fmt='%d')
        
        # GO term counts
        go_counts = [comparison_data['species1_unique_go_terms'], 
                    comparison_data['species2_unique_go_terms']]
        
        bars2 = ax2.bar(species, go_counts, color=[self.colors['accent'], self.colors['success']])
        ax2.set_ylabel('Number of Unique GO Terms')
        ax2.set_title('GO Term Counts by Species')
        ax2.bar_label(bars2, fmt='%d')
        
        # Venn diagram style comparison
        common = comparison_data['common_go_terms']
        unique1 = comparison_data['species1_unique_go_terms'] - common
        unique2 = comparison_data['species2_unique_go_terms'] - common
        
        venn_data = [unique1, unique2, common]
        venn_labels = ['Species 1\nUnique', 'Species 2\nUnique', 'Common']
        
        bars3 = ax3.bar(venn_labels, venn_data, 
                       color=[self.colors['primary'], self.colors['secondary'], self.colors['accent']])
        ax3.set_ylabel('Number of GO Terms')
        ax3.set_title('GO Term Overlap Analysis')
        ax3.bar_label(bars3, fmt='%d')
        
        # Coverage metrics
        coverage_metrics = ['GO Term\nOverlap']
        coverage_values = [comparison_data['go_term_overlap']]
        
        bars4 = ax4.bar(coverage_metrics, coverage_values, color=self.colors['neutral'])
        ax4.set_ylabel('Coverage Ratio')
        ax4.set_title('Coverage Metrics')
        ax4.set_ylim(0, 1)
        ax4.bar_label(bars4, fmt='%.3f')
        
        plt.suptitle('Species Comparison Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Species comparison plot saved to {save_path}")
        
        plt.show()
    
    def plot_alignment_quality_heatmap(self, detailed_scores: pd.DataFrame,
                                     save_path: Optional[Path] = None) -> None:
        """
        Create a heatmap showing alignment quality patterns.
        
        Args:
            detailed_scores: DataFrame with detailed pair scores
            save_path: Path to save the plot
        """
        # Create bins for GO term counts
        detailed_scores['go_terms_1_bin'] = pd.cut(detailed_scores['go_terms_1_count'], 
                                                 bins=5, labels=['Very Low', 'Low', 'Medium', 'High', 'Very High'])
        detailed_scores['go_terms_2_bin'] = pd.cut(detailed_scores['go_terms_2_count'], 
                                                 bins=5, labels=['Very Low', 'Low', 'Medium', 'High', 'Very High'])
        
        # Create pivot table for heatmap
        heatmap_data = detailed_scores.groupby(['go_terms_1_bin', 'go_terms_2_bin'])['similarity_score'].mean().unstack()
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='YlOrRd', 
                   cbar_kws={'label': 'Average Similarity Score'})
        plt.title('Alignment Quality Heatmap\nby GO Term Count Ranges', fontsize=14, fontweight='bold')
        plt.xlabel('Species 2 GO Terms')
        plt.ylabel('Species 1 GO Terms')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Quality heatmap saved to {save_path}")
        
        plt.show()
    
    def create_interactive_dashboard(self, results: Dict[str, any], 
                                   detailed_scores: pd.DataFrame,
                                   save_path: Optional[Path] = None) -> None:
        """
        Create an interactive Plotly dashboard.
        
        Args:
            results: Results dictionary from scoring
            detailed_scores: DataFrame with detailed scores
            save_path: Path to save the HTML dashboard
        """
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Score Distribution', 'Species Comparison', 
                          'Score vs GO Terms', 'Quality Distribution'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Score distribution histogram
        fig.add_trace(
            go.Histogram(x=detailed_scores['similarity_score'], 
                        nbinsx=30, name='Score Distribution',
                        marker_color=self.colors['primary']),
            row=1, col=1
        )
        
        # Species comparison bar chart
        comparison = results['species_comparison']
        species = ['Species 1', 'Species 2']
        protein_counts = [comparison['species1_total_proteins'], 
                         comparison['species2_total_proteins']]
        
        fig.add_trace(
            go.Bar(x=species, y=protein_counts, name='Protein Counts',
                  marker_color=[self.colors['primary'], self.colors['secondary']]),
            row=1, col=2
        )
        
        # Score vs GO terms scatter plot
        fig.add_trace(
            go.Scatter(x=detailed_scores['go_terms_1_count'], 
                      y=detailed_scores['similarity_score'],
                      mode='markers', name='Species 1',
                      marker=dict(color=self.colors['primary'], opacity=0.6)),
            row=2, col=1
        )
        
        fig.add_trace(
            go.Scatter(x=detailed_scores['go_terms_2_count'], 
                      y=detailed_scores['similarity_score'],
                      mode='markers', name='Species 2',
                      marker=dict(color=self.colors['secondary'], opacity=0.6)),
            row=2, col=1
        )
        
        # Quality distribution box plot
        fig.add_trace(
            go.Box(y=detailed_scores['similarity_score'], name='Quality Distribution',
                  marker_color=self.colors['accent']),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            title_text="Network Alignment Quality Dashboard",
            showlegend=True,
            height=800
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Similarity Score", row=1, col=1)
        fig.update_yaxes(title_text="Frequency", row=1, col=1)
        fig.update_xaxes(title_text="Species", row=1, col=2)
        fig.update_yaxes(title_text="Protein Count", row=1, col=2)
        fig.update_xaxes(title_text="GO Terms Count", row=2, col=1)
        fig.update_yaxes(title_text="Similarity Score", row=2, col=1)
        fig.update_xaxes(title_text="Quality Distribution", row=2, col=2)
        fig.update_yaxes(title_text="Similarity Score", row=2, col=2)
        
        if save_path:
            fig.write_html(save_path)
            logger.info(f"Interactive dashboard saved to {save_path}")
        
        fig.show()
    
    def plot_go_term_analysis(self, go_statistics: Dict[str, any],
                            save_path: Optional[Path] = None) -> None:
        """
        Plot GO term analysis and statistics.
        
        Args:
            go_statistics: GO statistics from analyzer
            save_path: Path to save the plot
        """
        if not go_statistics:
            logger.warning("No GO statistics provided")
            return
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Most common GO terms
        if 'most_common_go_terms' in go_statistics:
            terms, counts = zip(*go_statistics['most_common_go_terms'][:10])
            bars = ax1.barh(range(len(terms)), counts, color=self.colors['primary'])
            ax1.set_yticks(range(len(terms)))
            ax1.set_yticklabels(terms, fontsize=8)
            ax1.set_xlabel('Annotation Count')
            ax1.set_title('Most Common GO Terms')
            ax1.bar_label(bars, fmt='%d')
        
        # Annotation distribution
        if 'avg_annotations_per_protein' in go_statistics:
            metrics = ['Min', 'Avg', 'Max']
            values = [go_statistics['min_annotations_per_protein'],
                     go_statistics['avg_annotations_per_protein'],
                     go_statistics['max_annotations_per_protein']]
            
            bars = ax2.bar(metrics, values, color=[self.colors['success'], 
                                                 self.colors['accent'], 
                                                 self.colors['secondary']])
            ax2.set_ylabel('Annotations per Protein')
            ax2.set_title('Annotation Distribution')
            ax2.bar_label(bars, fmt='%.1f')
        
        # Coverage pie chart
        if 'total_proteins' in go_statistics and 'total_annotations' in go_statistics:
            annotated = go_statistics['total_annotations']
            total = go_statistics['total_proteins']
            unannotated = total - annotated
            
            ax3.pie([annotated, unannotated], 
                   labels=['Annotated', 'Unannotated'],
                   colors=[self.colors['primary'], self.colors['neutral']],
                   autopct='%1.1f%%')
            ax3.set_title('Protein Annotation Coverage')
        
        # GO term diversity
        if 'unique_go_terms' in go_statistics:
            ax4.text(0.5, 0.5, f"Unique GO Terms:\n{go_statistics['unique_go_terms']:,}",
                    ha='center', va='center', fontsize=20, fontweight='bold',
                    transform=ax4.transAxes)
            ax4.set_title('GO Term Diversity')
            ax4.axis('off')
        
        plt.suptitle('Gene Ontology Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"GO analysis plot saved to {save_path}")
        
        plt.show() 