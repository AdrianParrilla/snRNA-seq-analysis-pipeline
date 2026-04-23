#!/usr/bin/env python3
"""
Validate cellbender samplesheet format.
Checks required columns and file existence.
"""

import sys
import csv
import re
import argparse
from pathlib import Path
import pandas as pd
from typing import List, Tuple

def check_samplesheet(samplesheet_file: str):
    """
    Validate cellbendersamplesheet format.

    Args:
        samplesheet_file: Path to TSV samplesheet

    Returns:
        (is_valid, error_messages)
    """

    errors = []

    try:
        samplesheet = pd.read_csv(samplesheet_file, sep='\t')

        # Check required columns
        required_cols = ['sample', 'matrix_path']
        missing_cols = [col for col in required_cols if col not in samplesheet.columns]

        if missing_cols:
            errors.append(f"Missing required columns: {missing}")
            return False, errors

        # Validate each row
        for row_num, row in samplesheet.iterrows():
            sample_id = row['sample'].strip()
            matrix_path = row['matrix_path'].strip()
            learning_rate = row['learning_rate']
            expected_cells = row['expected_cells']
            total_droplets_included = row['total_droplets_included']

            # Validate sample_id
            if not sample_id:
                errors.append(f"Row {row_num}: sample_id is empty")
            elif not re.match(r'^[a-zA-Z0-9_\-\.]+$', sample_id):
                errors.append(f"Row {row_num}: Invalid sample_id format: {sample_id}")

            # Validate matrix_path
            if not matrix_path:
                errors.append(f"Row {row_num}: chemistry is empty")

    except Exception as e:
        errors.append(f"Unexpected error: {e}")
        return False, errors

    return len(errors) == 0, errors


def main(samplesheet):
    
    is_valid, errors = check_samplesheet(samplesheet)

    if is_valid:
        print(f"✓ Samplesheet '{args.samplesheet}' is valid")
        sys.exit(0)
    else:
        print(f"✗ Samplesheet '{args.samplesheet}' has errors:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Validate cellbender samplesheet')
    parser.add_argument('samplesheet', help='Samplesheet TSV file')

    args = parser.parse_args()

    main(args.samplesheet)