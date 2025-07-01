"""
Similarity Metrics Module

Provides various similarity measures for comparing proteins based on GO annotations.
"""

import logging
from typing import Set, Dict, List, Optional, Tuple
import numpy as np
from scipy.spatial.distance import jaccard
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict

logger = logging.getLogger(__name__)


class JaccardSimilarity:
    """Jaccard similarity calculator for GO term sets."""
    
    @staticmethod
    def calculate(go_terms_1: Set[str], go_terms_2: Set[str]) -> float:
        """
        Calculate Jaccard similarity between two sets of GO terms.
        
        Args:
            go_terms_1: First set of GO terms
            go_terms_2: Second set of GO terms
            
        Returns:
            Jaccard similarity score (0-1)
        """
        if not go_terms_1 and not go_terms_2:
            return 1.0  # Both empty sets are considered identical
        
        if not go_terms_1 or not go_terms_2:
            return 0.0  # One empty set means no similarity
        
        intersection = len(go_terms_1.intersection(go_terms_2))
        union = len(go_terms_1.union(go_terms_2))
        
        return intersection / union if union > 0 else 0.0


class SemanticSimilarity:
    """Semantic similarity calculator using GO term hierarchies."""
    
    def __init__(self, go_hierarchy: Optional[Dict[str, Set[str]]] = None):
        """
        Initialize semantic similarity calculator.
        
        Args:
            go_hierarchy: Dictionary mapping GO terms to their parent terms
        """
        self.go_hierarchy = go_hierarchy or {}
        
    def calculate_resnik_similarity(self, go_terms_1: Set[str], 
                                  go_terms_2: Set[str]) -> float:
        """
        Calculate Resnik semantic similarity.
        
        Args:
            go_terms_1: First set of GO terms
            go_terms_2: Second set of GO terms
            
        Returns:
            Resnik similarity score
        """
        if not go_terms_1 or not go_terms_2:
            return 0.0
            
        # For now, implement a simplified version
        # In a full implementation, you would need GO term information content
        max_similarity = 0.0
        
        for term1 in go_terms_1:
            for term2 in go_terms_2:
                similarity = self._calculate_term_similarity(term1, term2)
                max_similarity = max(max_similarity, similarity)
        
        return max_similarity
    
    def _calculate_term_similarity(self, term1: str, term2: str) -> float:
        """
        Calculate similarity between two individual GO terms.
        
        Args:
            term1: First GO term
            term2: Second GO term
            
        Returns:
            Similarity score
        """
        # Simplified implementation - in practice you'd use GO term hierarchies
        if term1 == term2:
            return 1.0
        
        # Check if terms are related in hierarchy
        if self.go_hierarchy:
            if term1 in self.go_hierarchy.get(term2, set()) or term2 in self.go_hierarchy.get(term1, set()):
                return 0.5
        
        return 0.0


class CosineSimilarity:
    """Cosine similarity calculator for GO term vectors."""
    
    @staticmethod
    def calculate(go_terms_1: Set[str], go_terms_2: Set[str], 
                 all_terms: Optional[Set[str]] = None) -> float:
        """
        Calculate cosine similarity between GO term sets.
        
        Args:
            go_terms_1: First set of GO terms
            go_terms_2: Second set of GO terms
            all_terms: Set of all possible GO terms (for vectorization)
            
        Returns:
            Cosine similarity score (0-1)
        """
        if not all_terms:
            all_terms = go_terms_1.union(go_terms_2)
        
        if not all_terms:
            return 1.0  # Both empty
        
        # Create binary vectors
        vector1 = np.array([1 if term in go_terms_1 else 0 for term in all_terms])
        vector2 = np.array([1 if term in go_terms_2 else 0 for term in all_terms])
        
        # Calculate cosine similarity
        dot_product = np.dot(vector1, vector2)
        norm1 = np.linalg.norm(vector1)
        norm2 = np.linalg.norm(vector2)
        
        if norm1 == 0 or norm2 == 0:
            return 0.0
        
        return dot_product / (norm1 * norm2)


class AlignmentMetrics:
    """Comprehensive metrics for network alignment evaluation."""
    
    def __init__(self):
        self.metrics = {}
        
    def calculate_alignment_score(self, alignment_pairs: List[Tuple[str, str]], 
                                go_annotations_1: Dict[str, Set[str]], 
                                go_annotations_2: Dict[str, Set[str]], 
                                metric: str = "jaccard") -> Dict[str, float]:
        """
        Calculate comprehensive alignment scores.
        
        Args:
            alignment_pairs: List of (protein1, protein2) pairs
            go_annotations_1: GO annotations for first species
            go_annotations_2: GO annotations for second species
            metric: Similarity metric to use ("jaccard", "cosine", "resnik")
            
        Returns:
            Dictionary with various alignment metrics
        """
        scores = []
        unmappable_1 = 0
        unmappable_2 = 0
        
        for protein1, protein2 in alignment_pairs:
            go_terms_1 = go_annotations_1.get(protein1, set())
            go_terms_2 = go_annotations_2.get(protein2, set())
            
            if not go_terms_1:
                unmappable_1 += 1
            if not go_terms_2:
                unmappable_2 += 1
            
            if go_terms_1 and go_terms_2:
                if metric == "jaccard":
                    score = JaccardSimilarity.calculate(go_terms_1, go_terms_2)
                elif metric == "cosine":
                    score = CosineSimilarity.calculate(go_terms_1, go_terms_2)
                elif metric == "resnik":
                    semantic_sim = SemanticSimilarity()
                    score = semantic_sim.calculate_resnik_similarity(go_terms_1, go_terms_2)
                else:
                    raise ValueError(f"Unknown metric: {metric}")
                
                scores.append(score)
        
        if not scores:
            return {
                "mean_score": 0.0,
                "median_score": 0.0,
                "std_score": 0.0,
                "total_pairs": len(alignment_pairs),
                "scored_pairs": 0,
                "unmappable_1": unmappable_1,
                "unmappable_2": unmappable_2
            }
        
        scores_array = np.array(scores)
        
        return {
            "mean_score": float(np.mean(scores_array)),
            "median_score": float(np.median(scores_array)),
            "std_score": float(np.std(scores_array)),
            "min_score": float(np.min(scores_array)),
            "max_score": float(np.max(scores_array)),
            "total_pairs": len(alignment_pairs),
            "scored_pairs": len(scores),
            "unmappable_1": unmappable_1,
            "unmappable_2": unmappable_2,
            "coverage": len(scores) / len(alignment_pairs) if alignment_pairs else 0.0
        }
    
    def calculate_species_comparison(self, go_annotations_1: Dict[str, Set[str]], 
                                   go_annotations_2: Dict[str, Set[str]]) -> Dict[str, any]:
        """
        Calculate comparison statistics between two species.
        
        Args:
            go_annotations_1: GO annotations for first species
            go_annotations_2: GO annotations for second species
            
        Returns:
            Dictionary with comparison statistics
        """
        all_terms_1 = set().union(*go_annotations_1.values())
        all_terms_2 = set().union(*go_annotations_2.values())
        
        common_terms = all_terms_1.intersection(all_terms_2)
        
        return {
            "species1_total_proteins": len(go_annotations_1),
            "species2_total_proteins": len(go_annotations_2),
            "species1_unique_go_terms": len(all_terms_1),
            "species2_unique_go_terms": len(all_terms_2),
            "common_go_terms": len(common_terms),
            "go_term_overlap": len(common_terms) / len(all_terms_1.union(all_terms_2)) if all_terms_1.union(all_terms_2) else 0.0
        } 