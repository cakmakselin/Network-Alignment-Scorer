"""
Mapping Parser Module

Handles parsing of protein ID mapping files with support for multiple ID types.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set
import pandas as pd
from pydantic import BaseModel, validator

logger = logging.getLogger(__name__)


class MappingConfig(BaseModel):
    """Configuration for mapping file parsing."""
    
    delimiter: str = "\t"
    skip_header: bool = True
    ensembl_col: str = "Ensembl_ID"
    
    @validator("delimiter")
    def validate_delimiter(cls, v):
        if v not in ["\t", ",", ";", " "]:
            raise ValueError("Delimiter must be tab, comma, semicolon, or space")
        return v


class MappingParser:
    """
    Parser for protein ID mapping files.
    
    Supports various ID types and provides efficient lookup capabilities.
    """
    
    def __init__(self, config: Optional[MappingConfig] = None):
        self.config = config or MappingConfig()
        self.mappings: Dict[str, Dict[str, str]] = {}
        self.reverse_mappings: Dict[str, Dict[str, Set[str]]] = {}
        
    def parse_mapping_file(self, file_path: Path) -> Dict[str, Dict[str, str]]:
        """
        Parse a mapping file and return organized mappings.
        
        Args:
            file_path: Path to the mapping file
            
        Returns:
            Dictionary mapping ID types to their mappings
        """
        if not file_path.exists():
            raise FileNotFoundError(f"Mapping file not found: {file_path}")
            
        logger.info(f"Parsing mapping file: {file_path}")
        
        try:
            # Read the file
            df = pd.read_csv(
                file_path, 
                delimiter=self.config.delimiter,
                skiprows=1 if self.config.skip_header else 0
            )
            
            # Clean column names
            df.columns = df.columns.str.strip()
            
            # Ensure Ensembl ID column exists
            if self.config.ensembl_col not in df.columns:
                raise ValueError(f"Ensembl ID column '{self.config.ensembl_col}' not found")
            
            # Create mappings for each ID type
            mappings = {}
            for col in df.columns:
                if col != self.config.ensembl_col:
                    # Create mapping from this ID type to Ensembl
                    id_mapping = {}
                    for _, row in df.iterrows():
                        ensembl_id = row[self.config.ensembl_col]
                        other_id = row[col]
                        
                        if pd.notna(other_id) and pd.notna(ensembl_id):
                            id_mapping[str(other_id).strip()] = str(ensembl_id).strip()
                    
                    mappings[col] = id_mapping
                    
                    # Create reverse mapping for efficient lookup
                    reverse_mapping = {}
                    for other_id, ensembl_id in id_mapping.items():
                        if ensembl_id not in reverse_mapping:
                            reverse_mapping[ensembl_id] = set()
                        reverse_mapping[ensembl_id].add(other_id)
                    
                    self.reverse_mappings[col] = reverse_mapping
            
            self.mappings = mappings
            
            # Log statistics
            total_mappings = sum(len(mapping) for mapping in mappings.values())
            logger.info(f"Parsed {len(mappings)} ID types with {total_mappings} total mappings")
            
            return mappings
            
        except Exception as e:
            logger.error(f"Error parsing mapping file: {e}")
            raise
    
    def get_ensembl_id(self, id_value: str, id_type: str) -> Optional[str]:
        """
        Get Ensembl ID for a given ID value and type.
        
        Args:
            id_value: The ID value to look up
            id_type: The type of ID (column name)
            
        Returns:
            Ensembl ID if found, None otherwise
        """
        if id_type not in self.mappings:
            logger.warning(f"Unknown ID type: {id_type}")
            return None
            
        return self.mappings[id_type].get(id_value)
    
    def get_all_ids(self, ensembl_id: str, id_type: str) -> Set[str]:
        """
        Get all IDs of a specific type for a given Ensembl ID.
        
        Args:
            ensembl_id: The Ensembl ID
            id_type: The type of ID to retrieve
            
        Returns:
            Set of IDs of the specified type
        """
        if id_type not in self.reverse_mappings:
            return set()
            
        return self.reverse_mappings[id_type].get(ensembl_id, set())
    
    def get_mapping_stats(self) -> Dict[str, int]:
        """
        Get statistics about the parsed mappings.
        
        Returns:
            Dictionary with mapping statistics
        """
        stats = {}
        for id_type, mapping in self.mappings.items():
            stats[id_type] = len(mapping)
        return stats 