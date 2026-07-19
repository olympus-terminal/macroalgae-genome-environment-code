#!/usr/bin/env python3
"""
Extract AlphaEarth embeddings from Google Earth Engine for clean PFAM dataset.
Uses ONLY trusted data: pfam_counts_with_metadata_20251019.csv

Created: 2025-10-19
"""

import ee
import pandas as pd
import numpy as np
from datetime import datetime
import time

# Configuration - ONLY CLEAN DATA
INPUT_FILE = "pfam_counts_with_metadata_20251019.csv"
OUTPUT_FILE = f"alphaearth_embeddings_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"

# PREVENT ACCIDENTAL USE OF CORRUPTED DATA
FORBIDDEN_FILES = ["Table_S4_meta-and_pfams.csv"]
import os
for forbidden in FORBIDDEN_FILES:
    if os.path.exists(forbidden):
        raise ValueError(f"CORRUPTED DATA DETECTED: {forbidden} must be deleted before running!")

def initialize_gee():
    """Initialize Google Earth Engine."""
    print("Initializing Google Earth Engine...")
    try:
        ee.Initialize()
        print("✓ GEE initialized successfully")
    except Exception as e:
        print(f"✗ GEE initialization failed: {e}")
        print("\nPlease run: earthengine authenticate")
        raise

def load_clean_data():
    """Load clean PFAM data and prepare coordinates."""
    print(f"\nLoading clean data from {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE)

    # Filter for genomes with valid coordinates
    valid_coords = df['DD latitude'].notna() & df['DD longitude'].notna()
    df = df[valid_coords].copy()

    # Add ID number for tracking
    df['ID number'] = range(1, len(df) + 1)

    print(f"  Total genomes: {len(df)}")
    print(f"  Genomes with valid coordinates: {len(df)}")

    return df

def extract_embeddings(df):
    """Extract AlphaEarth embeddings for all genomes using batch method."""
    print("\n" + "=" * 70)
    print("Extracting AlphaEarth embeddings from Google Earth Engine...")
    print("=" * 70)

    try:
        # Load AlphaEarth collection and create mosaic
        ALPHAEARTH_COLLECTION = "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL"
        print(f"Loading collection: {ALPHAEARTH_COLLECTION}")
        alphaearth = ee.ImageCollection(ALPHAEARTH_COLLECTION).mosaic()
        print("✓ Collection loaded and mosaicked")

        # Create feature collection from coordinates
        features = []
        for idx, row in df.iterrows():
            point = ee.Geometry.Point([row['DD longitude'], row['DD latitude']])
            features.append(ee.Feature(point, {
                'id': int(row['ID number']),
                'genome': row['Genome'],
                'species': row['Species']
            }))

        feature_collection = ee.FeatureCollection(features)
        print(f"✓ Created FeatureCollection with {len(df)} points")

        # Extract embeddings using batch method
        print("  Extracting embeddings (this may take 2-3 minutes)...")
        sampled = alphaearth.reduceRegions(
            collection=feature_collection,
            reducer=ee.Reducer.mean(),
            scale=10  # 10 meter scale
        )

        # Fetch results
        results = sampled.getInfo()
        print("✓ Retrieved results from GEE")

        # Process results
        embedding_data = []
        failed_ids = []

        for feature in results['features']:
            props = feature['properties']
            id_num = props['id']

            # Get corresponding genome info
            genome_row = df[df['ID number'] == id_num].iloc[0]

            row_data = {
                'Genome': genome_row['Genome'],
                'Species': genome_row['Species'],
                'ID number': id_num,
                'DD latitude': genome_row['DD latitude'],
                'DD longitude': genome_row['DD longitude']
            }

            # Extract embedding dimensions (A00-A63)
            has_data = False
            for i in range(64):
                band_name = f'A{i:02d}'
                value = props.get(band_name, None)
                if value is not None and not (isinstance(value, float) and np.isnan(value)):
                    has_data = True
                row_data[band_name] = value if value is not None else np.nan

            if has_data:
                embedding_data.append(row_data)
            else:
                failed_ids.append(id_num)

        # Convert to DataFrame
        if embedding_data:
            embedding_df = pd.DataFrame(embedding_data)

            # Ensure columns are in order: metadata first, then A00-A63
            metadata_cols = ['Genome', 'Species', 'ID number', 'DD latitude', 'DD longitude']
            embedding_cols = [f'A{i:02d}' for i in range(64)]
            embedding_df = embedding_df[metadata_cols + embedding_cols]

            print("\n" + "=" * 70)
            print("SUMMARY")
            print("=" * 70)
            print(f"  Total attempts: {len(df)}")
            print(f"  Successful: {len(embedding_data)}")
            print(f"  Failed: {len(failed_ids)}")
            print(f"  Success rate: {100 * len(embedding_data) / len(df):.1f}%")
            print(f"  Embedding dimensions: {len(embedding_cols)}")

            if failed_ids:
                print(f"\nFailed IDs: {failed_ids}")
                failed_genomes = df[df['ID number'].isin(failed_ids)]
                for _, row in failed_genomes.iterrows():
                    print(f"  - {row['Species']} ({row['DD latitude']:.2f}, {row['DD longitude']:.2f})")

            return embedding_df
        else:
            raise ValueError("No embeddings were successfully extracted!")

    except Exception as e:
        print(f"\n✗ Extraction failed: {e}")
        raise

def main():
    """Main execution function."""
    print("=" * 70)
    print("AlphaEarth Embedding Extraction")
    print("=" * 70)
    print(f"Input: {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")

    # Initialize GEE
    initialize_gee()

    # Load data
    df = load_clean_data()

    # Extract embeddings
    embedding_df = extract_embeddings(df)

    # Save results
    print(f"\nSaving embeddings to {OUTPUT_FILE}...")
    embedding_df.to_csv(OUTPUT_FILE, index=False)
    print(f"✓ Saved {len(embedding_df)} embeddings")

    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print(f"\nOutput file: {OUTPUT_FILE}")
    print(f"Shape: {embedding_df.shape}")

if __name__ == "__main__":
    main()
