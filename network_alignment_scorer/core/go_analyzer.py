"""
Gene Ontology Analyzer Module

Handles parsing and analysis of Gene Ontology annotation files.
"""

import logging
from pathlib import Path
from typing import Dict, Set, Optional, List, Tuple
import pandas as pd
from collections import defaultdict
from pydantic import BaseModel

logger = logging.getLogger(__name__)


class GOConfig(BaseModel):
    """Configuration for GO file parsing."""
    
    comment_char: str = "!"
    go_prefix: str = "GO:"
    id_column: int = 1
    go_column: Optional[int] = None  # Auto-detect if None
    
    
class GOAnalyzer:
    """
    Analyzer for Gene Ontology annotation files.
    
    Provides efficient parsing and lookup of GO terms for proteins.
    """
    
    def __init__(self, config: Optional[GOConfig] = None):
        self.config = config or GOConfig()
        self.go_annotations: Dict[str, Set[str]] = {}
        self.protein_stats: Dict[str, int] = defaultdict(int)
        self.go_term_stats: Dict[str, int] = defaultdict(int)
        
    def parse_go_file(self, file_path: Path, mapping_parser) -> Dict[str, Set[str]]:
        """
        Parse a GO annotation file and return protein-to-GO mappings.
        
        Args:
            file_path: Path to the GO annotation file
            mapping_parser: MappingParser instance for ID conversion
            
        Returns:
            Dictionary mapping Ensembl IDs to sets of GO terms
        """
        if not file_path.exists():
            raise FileNotFoundError(f"GO file not found: {file_path}")
            
        logger.info(f"Parsing GO file: {file_path}")
        
        go_annotations = defaultdict(set)
        line_count = 0
        mapped_count = 0
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    line_count += 1
                    
                    # Skip comments and empty lines
                    if not line or line.startswith(self.config.comment_char):
                        continue
                    
                    # Parse the line
                    parts = line.split()
                    if len(parts) < 2:
                        logger.warning(f"Skipping malformed line {line_count}: {line}")
                        continue
                    
                    # Extract protein ID and GO terms
                    protein_id = parts[self.config.id_column]
                    
                    # Find GO terms in the line
                    go_terms = []
                    for part in parts:
                        if part.startswith(self.config.go_prefix):
                            go_term = part[len(self.config.go_prefix):]  # Remove prefix
                            go_terms.append(go_term)
                    
                    if not go_terms:
                        logger.warning(f"No GO terms found in line {line_count}: {line}")
                        continue
                    
                    # Try to map protein ID to Ensembl ID
                    ensembl_id = None
                    for id_type in mapping_parser.mappings.keys():
                        ensembl_id = mapping_parser.get_ensembl_id(protein_id, id_type)
                        if ensembl_id:
                            break
                    
                    if ensembl_id:
                        # Add GO terms to the protein
                        for go_term in go_terms:
                            go_annotations[ensembl_id].add(go_term)
                            self.go_term_stats[go_term] += 1
                        
                        self.protein_stats[ensembl_id] += len(go_terms)
                        mapped_count += 1
                    else:
                        logger.debug(f"Could not map protein ID: {protein_id}")
            
            self.go_annotations = dict(go_annotations)
            
            # Log statistics
            total_go_terms = sum(len(terms) for terms in go_annotations.values())
            logger.info(f"Parsed {line_count} lines, mapped {mapped_count} proteins")
            logger.info(f"Total GO annotations: {total_go_terms}")
            logger.info(f"Unique GO terms: {len(self.go_term_stats)}")
            
            return self.go_annotations
            
        except Exception as e:
            logger.error(f"Error parsing GO file: {e}")
            raise
    
    def get_go_terms(self, protein_id: str) -> Set[str]:
        """
        Get GO terms for a specific protein.
        
        Args:
            protein_id: Ensembl ID of the protein
            
        Returns:
            Set of GO terms associated with the protein
        """
        return self.go_annotations.get(protein_id, set())
    
    def get_proteins_with_go_term(self, go_term: str) -> Set[str]:
        """
        Get all proteins annotated with a specific GO term.
        
        Args:
            go_term: The GO term to search for
            
        Returns:
            Set of protein IDs with this GO term
        """
        proteins = set()
        for protein_id, go_terms in self.go_annotations.items():
            if go_term in go_terms:
                proteins.add(protein_id)
        return proteins
    
    def get_annotation_statistics(self) -> Dict[str, any]:
        """
        Get comprehensive statistics about the GO annotations.
        
        Returns:
            Dictionary with various statistics
        """
        if not self.go_annotations:
            return {}
        
        protein_go_counts = [len(terms) for terms in self.go_annotations.values()]
        
        stats = {
            "total_proteins": len(self.go_annotations),
            "total_annotations": sum(protein_go_counts),
            "unique_go_terms": len(self.go_term_stats),
            "avg_annotations_per_protein": sum(protein_go_counts) / len(protein_go_counts),
            "min_annotations_per_protein": min(protein_go_counts),
            "max_annotations_per_protein": max(protein_go_counts),
            "most_common_go_terms": sorted(
                self.go_term_stats.items(), 
                key=lambda x: x[1], 
                reverse=True
            )[:10]
        }
        
        return stats
    
    def filter_by_go_terms(self, required_terms: Set[str], 
                          min_terms: int = 1) -> Dict[str, Set[str]]:
        """
        Filter proteins to only those with at least min_terms from required_terms.
        
        Args:
            required_terms: Set of GO terms to check for
            min_terms: Minimum number of required terms (default: 1)
            
        Returns:
            Filtered dictionary of protein to GO term mappings
        """
        filtered = {}
        
        for protein_id, go_terms in self.go_annotations.items():
            common_terms = go_terms.intersection(required_terms)
            if len(common_terms) >= min_terms:
                filtered[protein_id] = go_terms
        
        return filtered 