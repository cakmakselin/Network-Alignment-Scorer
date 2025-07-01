"""
Network Alignment Scorer

Main module for scoring protein-protein interaction network alignments.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
import numpy as np
from collections import defaultdict

logger = logging.getLogger(__name__)


class NetworkAlignmentScorer:
    """
    Main class for scoring network alignments using GO annotations.
    
    Provides evaluation of protein-protein interaction network alignments
    across different species using Gene Ontology annotations.
    """
    
    def __init__(self):
        """Initialize the network alignment scorer."""
        self.mappings: Dict[str, Dict[str, str]] = {}
        self.go_annotations: Dict[str, Dict[str, Set[str]]] = {}
        self.alignments: List[Tuple[str, str]] = []
        
    def get_mapping(self, map_file: Path) -> List[Dict[str, str]]:
        """
        Parse mapping file and return list of dictionaries.
        This matches the original assignment logic exactly.
        
        Args:
            map_file: Path to mapping file
            
        Returns:
            List of dictionaries mapping different ID types to Ensembl IDs
        """
        logger.info(f"Parsing mapping file: {map_file}")
        
        if not map_file.exists():
            raise FileNotFoundError(f"Mapping file not found: {map_file}")
        
        mapping_list = []
        
        try:
            with open(map_file, 'r') as f:
                # Skip the header on the first line
                header = f.readline()
                
                # Create dictionaries (# of non-Ensembl values)
                n = len(header.strip().split()) - 1
                while n > 0:
                    mapping_list.append({})
                    n = n - 1
                
                # Fill dictionaries
                for line in f:
                    if line == header:
                        continue
                    else:
                        # Strip line of whitespace, split elements into list
                        x = line.strip().split()
                        
                        # Assign first item (Ensembl_ID) as value
                        value = x[0]
                        
                        # Fill n-th dictionary with n+1th item as key
                        for key in list(range(1, len(x))):
                            if len(list(range(len(x)))) == len(header.strip().split()):
                                # Select correct (n-th) dictionary position in mapping_list
                                position = key - 1
                                # Store key and value in n-th dictionary in mapping_list
                                mapping_list[position].update({x[key]: value})
                            elif len(list(range(len(x)))) == 1:
                                break
            
            logger.info(f"Parsed {len(mapping_list)} ID types")
            return mapping_list
            
        except Exception as e:
            logger.error(f"Error parsing mapping file: {e}")
            raise
    
    def get_go_terms(self, mapping_list: List[Dict[str, str]], go_file: Path) -> Dict[str, Set[str]]:
        """
        Parse GO file and return protein-to-GO mappings.
        This matches the original assignment logic exactly.
        
        Args:
            mapping_list: List of mapping dictionaries
            go_file: Path to GO annotation file
            
        Returns:
            Dictionary mapping Ensembl IDs to sets of GO terms
        """
        logger.info(f"Parsing GO file: {go_file}")
        
        if not go_file.exists():
            raise FileNotFoundError(f"GO file not found: {go_file}")
        
        go_dict = dict()
        
        try:
            with open(go_file, 'r') as f:
                for line in f:
                    if line.startswith("!"):
                        continue  # Skip comments
                    else:
                        # Strip line of go-file from whitespace and split into list
                        x = line.strip().split()
                        
                        # Use ID in go file as search word
                        query = x[1]
                        
                        # Find GO term and save as value (original bug: only gets last GO term)
                        for item in x:
                            if item.startswith("GO:"):
                                value = item.strip("GO:")
                            else:
                                continue
                        
                        # Search for query in every dictionary in mapping list
                        for d in mapping_list:
                            if query in d.keys():
                                # If found, save Ensembl_ID as key
                                key = mapping_list[mapping_list.index(d)][query]
                                
                                # Update dictionary with key:value
                                if key in go_dict:
                                    go_dict[key].add(value)
                                else:
                                    go_dict[key] = {value}
                            else:
                                continue
            
            logger.info(f"Parsed GO annotations for {len(go_dict)} proteins")
            return go_dict
            
        except Exception as e:
            logger.error(f"Error parsing GO file: {e}")
            raise
    
    def compute_score(self, alignment_file: Path, go_one_dict: Dict[str, Set[str]], 
                     go_two_dict: Dict[str, Set[str]]) -> Tuple[int, int, float]:
        """
        Compute alignment score using Jaccard similarity.
        This matches the original assignment logic exactly.
        
        Args:
            alignment_file: Path to alignment file
            go_one_dict: GO annotations for first species
            go_two_dict: GO annotations for second species
            
        Returns:
            Tuple of (unmappable_1, unmappable_2, total_score)
        """
        logger.info(f"Computing alignment score: {alignment_file}")
        
        if not alignment_file.exists():
            raise FileNotFoundError(f"Alignment file not found: {alignment_file}")
        
        unmappable_one = 0
        unmappable_two = 0
        score = 0.0
        
        try:
            with open(alignment_file, 'r') as f:
                for line in f:
                    # Strip line from whitespace and split into list
                    x = line.strip().split()
                    
                    # Find go-terms for first protein from go_dictionary
                    protein1 = x[0]
                    
                    if protein1 in go_one_dict:
                        protein1_go = go_one_dict[protein1]
                    elif protein1 in go_two_dict:
                        protein1_go = go_two_dict[protein1]
                    else:
                        unmappable_one = unmappable_one + 1
                        continue  # Skip this line if protein1 is unmappable
                    
                    # Find go-terms for second protein from go_dictionary
                    protein2 = x[1]
                    
                    if protein2 in go_one_dict:
                        protein2_go = go_one_dict[protein2]
                    elif protein2 in go_two_dict:
                        protein2_go = go_two_dict[protein2]
                    else:
                        unmappable_two = unmappable_two + 1
                        continue  # Skip this line if protein2 is unmappable
                    
                    # Compute jaccard score for 2 proteins
                    line_score = len(protein1_go.intersection(protein2_go)) / len(protein1_go.union(protein2_go))
                    
                    # Add line score to global score for file
                    score = score + line_score
            
            logger.info(f"Computed score for alignment pairs")
            return unmappable_one, unmappable_two, score
            
        except Exception as e:
            logger.error(f"Error computing score: {e}")
            raise
    
    def score_alignment(self, alignment_file: Path, go_file_1: Path, go_file_2: Path,
                       mapping_file_1: Path, mapping_file_2: Path) -> Dict[str, any]:
        """
        Score a network alignment using GO annotations.
        This matches the original assignment logic exactly.
        
        Args:
            alignment_file: Path to alignment file (.sif format)
            go_file_1: Path to GO annotations for first species
            go_file_2: Path to GO annotations for second species
            mapping_file_1: Path to mapping file for first species
            mapping_file_2: Path to mapping file for second species
            
        Returns:
            Dictionary with scoring results
        """
        logger.info("Starting alignment scoring...")
        
        # Parse mapping files
        mapping_list_1 = self.get_mapping(mapping_file_1)
        mapping_list_2 = self.get_mapping(mapping_file_2)
        
        # Parse GO files
        go_dict_1 = self.get_go_terms(mapping_list_1, go_file_1)
        go_dict_2 = self.get_go_terms(mapping_list_2, go_file_2)
        
        # Compute score
        unmappable_1, unmappable_2, total_score = self.compute_score(
            alignment_file, go_dict_1, go_dict_2
        )
        
        # Calculate statistics (matching original logic)
        total_pairs = sum(1 for line in open(alignment_file) 
                         if line.strip() and not line.startswith('#'))
        scored_pairs = total_pairs - unmappable_1 - unmappable_2
        mean_score = total_score / scored_pairs if scored_pairs > 0 else 0.0
        
        results = {
            "total_pairs": total_pairs,
            "scored_pairs": scored_pairs,
            "unmappable_1": unmappable_1,
            "unmappable_2": unmappable_2,
            "total_score": total_score,
            "mean_score": mean_score,
            "coverage": scored_pairs / total_pairs if total_pairs > 0 else 0.0
        }
        
        logger.info(f"Alignment scoring completed. Mean score: {mean_score:.4f}")
        return results
    
    def get_quality_report(self, results: Dict[str, any]) -> str:
        """
        Generate a quality report for the alignment.
        
        Args:
            results: Results from score_alignment()
            
        Returns:
            Formatted quality report string
        """
        report = f"""
Network Alignment Quality Report
================================

Alignment Statistics:
- Total alignment pairs: {results['total_pairs']:,}
- Successfully scored pairs: {results['scored_pairs']:,}
- Coverage: {results['coverage']:.2%}
- Unmappable proteins (species1): {results['unmappable_1']:,}
- Unmappable proteins (species2): {results['unmappable_2']:,}

Similarity Scores:
- Total score: {results['total_score']:.4f}
- Mean similarity: {results['mean_score']:.4f}

Quality Assessment:
"""
        
        # Add quality assessment
        if results['mean_score'] >= 0.7:
            report += "- Overall Quality: EXCELLENT\n"
        elif results['mean_score'] >= 0.5:
            report += "- Overall Quality: GOOD\n"
        elif results['mean_score'] >= 0.3:
            report += "- Overall Quality: FAIR\n"
        else:
            report += "- Overall Quality: POOR\n"
            
        if results['coverage'] >= 0.8:
            report += "- Coverage: EXCELLENT\n"
        elif results['coverage'] >= 0.6:
            report += "- Coverage: GOOD\n"
        else:
            report += "- Coverage: NEEDS IMPROVEMENT\n"
            
        return report 