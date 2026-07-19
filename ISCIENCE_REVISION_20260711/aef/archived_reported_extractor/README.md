# Archived reported AEF extractor

`extract_aef_embeddings_20251019.py` is the actual historical script associated with the archived AEF embedding table used by the reported analyses.

- Workspace source: `DATA_S2_25DEC2025/extract_aef_embeddings_20251019.py`
- Source and released-file SHA-256: `3a1588c2d89b2a0f7e20afb333074b4e1a958b4c1dc5ebda6b89460277228d8a`
- Earth Engine collection: `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`
- Historical operation: `ee.ImageCollection(ALPHAEARTH_COLLECTION).mosaic()`
- Requested scale: 10 m

The script does not filter the collection to a calendar year before mosaicking. The archived output therefore has unfiltered annual-collection provenance, with overlap precedence determined by Earth Engine collection order. It must not be described as a year-specific extraction.

The later `../extract_exact_id_aef_embeddings_20260718_224936.py` is a separate prospective extractor that requires an explicit year and validates the exact 126-genome cohort. It was not used to create the archived reported AEF table.
