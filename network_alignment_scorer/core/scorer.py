"""
Network Alignment Scorer

Main module for scoring protein-protein interaction network alignments.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
from tqdm import tqdm

from .mapping_parser import MappingParser
from .go_analyzer import GOAnalyzer
from ..utils.metrics import AlignmentMetrics, JaccardSimilarity, CosineSimilarity

logger = logging.getLogger(__name__)


class NetworkAlignmentScorer:
    """
    Main class for scoring network alignments using GO annotations.
    
    Provides comprehensive evaluation of protein-protein interaction
    network alignments across different species.
    """
    
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the network alignment scorer.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config or {}
        self.mapping_parser = MappingParser()
        self.go_analyzer = GOAnalyzer()
        self.metrics = AlignmentMetrics()
        
        # Store parsed data
        self.mappings: Dict[str, Dict[str, Dict[str, str]]] = {}
        self.go_annotations: Dict[str, Dict[str, Set[str]]] = {}
        self.alignments: Dict[str, List[Tuple[str, str]]] = {}
        
    def load_data(self, alignment_file: Path, go_file_1: Path, go_file_2: Path,
                  mapping_file_1: Path, mapping_file_2: Path) -> None:
        """
        Load all required data files.
        
        Args:
            alignment_file: Path to alignment file (.sif format)
            go_file_1: Path to GO annotations for first species
            go_file_2: Path to GO annotations for second species
            mapping_file_1: Path to mapping file for first species
            mapping_file_2: Path to mapping file for second species
        """
        logger.info("Loading data files...")
        
        # Parse mapping files
        logger.info("Parsing mapping files...")
        self.mappings["species1"] = self.mapping_parser.parse_mapping_file(mapping_file_1)
        self.mappings["species2"] = self.mapping_parser.parse_mapping_file(mapping_file_2)
        
        # Parse GO files
        logger.info("Parsing GO annotation files...")
        self.go_annotations["species1"] = self.go_analyzer.parse_go_file(go_file_1, self.mapping_parser)
        self.go_analyzer = GOAnalyzer()  # Reset for second species
        self.go_annotations["species2"] = self.go_analyzer.parse_go_file(go_file_2, self.mapping_parser)
        
        # Parse alignment file
        logger.info("Parsing alignment file...")
        self.alignments["main"] = self._parse_alignment_file(alignment_file)
        
        logger.info("Data loading completed successfully")
        
    def _parse_alignment_file(self, alignment_file: Path) -> List[Tuple[str, str]]:
        """
        Parse alignment file in SIF format.
        
        Args:
            alignment_file: Path to alignment file
            
        Returns:
            List of (protein1, protein2) alignment pairs
        """
        if not alignment_file.exists():
            raise FileNotFoundError(f"Alignment file not found: {alignment_file}")
            
        alignments = []
        
        try:
            with open(alignment_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            protein1, protein2 = parts[0], parts[1]
                            alignments.append((protein1, protein2))
                            
        except Exception as e:
            logger.error(f"Error parsing alignment file: {e}")
            raise
            
        logger.info(f"Parsed {len(alignments)} alignment pairs")
        return alignments
    
    def score_alignment(self, metric: str = "jaccard") -> Dict[str, any]:
        """
        Score the loaded alignment using specified metric.
        
        Args:
            metric: Similarity metric to use ("jaccard", "cosine", "resnik")
            
        Returns:
            Dictionary with comprehensive scoring results
        """
        if not self.alignments.get("main"):
            raise ValueError("No alignment data loaded. Call load_data() first.")
            
        logger.info(f"Scoring alignment using {metric} similarity...")
        
        # Calculate alignment scores
        alignment_scores = self.metrics.calculate_alignment_score(
            self.alignments["main"],
            self.go_annotations["species1"],
            self.go_annotations["species2"],
            metric
        )
        
        # Calculate species comparison statistics
        species_comparison = self.metrics.calculate_species_comparison(
            self.go_annotations["species1"],
            self.go_annotations["species2"]
        )
        
        # Get GO annotation statistics
        go_stats_1 = self.go_analyzer.get_annotation_statistics()
        
        # Combine results
        results = {
            "alignment_scores": alignment_scores,
            "species_comparison": species_comparison,
            "go_statistics": {
                "species1": go_stats_1,
                "species2": self.go_analyzer.get_annotation_statistics()
            },
            "mapping_statistics": {
                "species1": self.mapping_parser.get_mapping_stats(),
                "species2": self.mapping_parser.get_mapping_stats()
            },
            "metadata": {
                "metric_used": metric,
                "total_alignments": len(self.alignments["main"]),
                "species1_proteins": len(self.go_annotations["species1"]),
                "species2_proteins": len(self.go_annotations["species2"])
            }
        }
        
        logger.info(f"Alignment scoring completed. Mean score: {alignment_scores['mean_score']:.4f}")
        
        return results
    
    def get_detailed_pair_scores(self, metric: str = "jaccard") -> pd.DataFrame:
        """
        Get detailed scores for each alignment pair.
        
        Args:
            metric: Similarity metric to use
            
        Returns:
            DataFrame with detailed pair-wise scores
        """
        if not self.alignments.get("main"):
            raise ValueError("No alignment data loaded. Call load_data() first.")
            
        results = []
        
        for protein1, protein2 in tqdm(self.alignments["main"], desc="Calculating pair scores"):
            go_terms_1 = self.go_annotations["species1"].get(protein1, set())
            go_terms_2 = self.go_annotations["species2"].get(protein2, set())
            
            if metric == "jaccard":
                score = JaccardSimilarity.calculate(go_terms_1, go_terms_2)
            elif metric == "cosine":
                score = CosineSimilarity.calculate(go_terms_1, go_terms_2)
            else:
                raise ValueError(f"Unknown metric: {metric}")
                
            results.append({
                "protein1": protein1,
                "protein2": protein2,
                "go_terms_1_count": len(go_terms_1),
                "go_terms_2_count": len(go_terms_2),
                "common_terms": len(go_terms_1.intersection(go_terms_2)),
                "union_terms": len(go_terms_1.union(go_terms_2)),
                "similarity_score": score,
                "has_annotation_1": bool(go_terms_1),
                "has_annotation_2": bool(go_terms_2)
            })
            
        return pd.DataFrame(results)
    
    def filter_high_quality_alignments(self, threshold: float = 0.5, 
                                     metric: str = "jaccard") -> List[Tuple[str, str]]:
        """
        Filter alignments to only high-quality pairs above threshold.
        
        Args:
            threshold: Minimum similarity score
            metric: Similarity metric to use
            
        Returns:
            List of high-quality alignment pairs
        """
        detailed_scores = self.get_detailed_pair_scores(metric)
        high_quality = detailed_scores[detailed_scores['similarity_score'] >= threshold]
        
        return list(zip(high_quality['protein1'], high_quality['protein2']))
    
    def export_results(self, output_dir: Path, results: Dict[str, any]) -> None:
        """
        Export results to various formats.
        
        Args:
            output_dir: Directory to save results
            results: Results dictionary from score_alignment()
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Export summary statistics
        summary_df = pd.DataFrame([results["alignment_scores"]])
        summary_df.to_csv(output_dir / "alignment_summary.csv", index=False)
        
        # Export detailed pair scores
        detailed_scores = self.get_detailed_pair_scores()
        detailed_scores.to_csv(output_dir / "detailed_pair_scores.csv", index=False)
        
        # Export species comparison
        species_df = pd.DataFrame([results["species_comparison"]])
        species_df.to_csv(output_dir / "species_comparison.csv", index=False)
        
        logger.info(f"Results exported to {output_dir}")
    
    def get_alignment_quality_report(self, results: Dict[str, any]) -> str:
        """
        Generate a human-readable quality report.
        
        Args:
            results: Results dictionary from score_alignment()
            
        Returns:
            Formatted quality report string
        """
        scores = results["alignment_scores"]
        comparison = results["species_comparison"]
        
        report = f"""
Network Alignment Quality Report
================================

Alignment Statistics:
- Total alignment pairs: {scores['total_pairs']:,}
- Successfully scored pairs: {scores['scored_pairs']:,}
- Coverage: {scores['coverage']:.2%}
- Unmappable proteins (species1): {scores['unmappable_1']:,}
- Unmappable proteins (species2): {scores['unmappable_2']:,}

Similarity Scores:
- Mean similarity: {scores['mean_score']:.4f}
- Median similarity: {scores['median_score']:.4f}
- Standard deviation: {scores['std_score']:.4f}
- Minimum score: {scores['min_score']:.4f}
- Maximum score: {scores['max_score']:.4f}

Species Comparison:
- Species1 proteins: {comparison['species1_total_proteins']:,}
- Species2 proteins: {comparison['species2_total_proteins']:,}
- Common GO terms: {comparison['common_go_terms']:,}
- GO term overlap: {comparison['go_term_overlap']:.2%}

Quality Assessment:
"""
        
        # Add quality assessment
        if scores['mean_score'] >= 0.7:
            report += "- Overall Quality: EXCELLENT\n"
        elif scores['mean_score'] >= 0.5:
            report += "- Overall Quality: GOOD\n"
        elif scores['mean_score'] >= 0.3:
            report += "- Overall Quality: FAIR\n"
        else:
            report += "- Overall Quality: POOR\n"
            
        if scores['coverage'] >= 0.8:
            report += "- Coverage: EXCELLENT\n"
        elif scores['coverage'] >= 0.6:
            report += "- Coverage: GOOD\n"
        else:
            report += "- Coverage: NEEDS IMPROVEMENT\n"
            
        return report 