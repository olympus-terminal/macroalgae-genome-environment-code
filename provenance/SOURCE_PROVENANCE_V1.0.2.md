# v1.0.2 source-script provenance

Paths in this table are relative to the private analysis-workspace root and reveal no host-specific directory. SHA-256 digests authenticate the exact workspace sources from which the public files were copied.

Most released files are byte-identical to their source. Seven supplemental builders received only a portable `#!/usr/bin/env python3` shebang in place of a host-specific interpreter path. The AEF–GEE alignment script received that source's computation unchanged and an affirmative docstring clarification: it is a descriptive same-site cross-representation alignment with distributed, nonunique axis meanings.

| Public release path | Workspace source path | Workspace source SHA-256 | Public-copy treatment |
|---|---|---|---|
| `ISCIENCE_REVISION_20260711/aef/archived_reported_extractor/extract_aef_embeddings_20251019.py` | `DATA_S2_25DEC2025/extract_aef_embeddings_20251019.py` | `3a1588c2d89b2a0f7e20afb333074b4e1a958b4c1dc5ebda6b89460277228d8a` | Byte-identical; explicitly archival and unfiltered |
| `ISCIENCE_REVISION_20260711/analysis_stats/correct_aef_structured_null_conditional_tail_20260719_205621.py` | Same relative path | `ef430c18d031d5e1ba2448505e60f7f8975d6991346d05c4d1913e8068a199ee` | Byte-identical |
| `ISCIENCE_REVISION_20260711/annotations/fetch_gee_robust_interpro_20260719_212508.py` | Same relative path | `b52c8d5f96feaf0e201678e3c5e10d30001126a049f43b4becf83175ac7cb8d3` | Byte-identical |
| `ISCIENCE_REVISION_20260711/figure5_gee_primary/build_figure5_primary_gee_sensitivity_20260719_210207.py` | Same relative path | `faa84795ff8f615efc4ef9552864aa3b0bafb4210a504b5f370e1a4a5de6bd48` | Byte-identical |
| `ISCIENCE_REVISION_20260711/gee_validation/run_aef_gee_site_alignment_20260719_204821.py` | Same relative path | `f514e8c6674ebea16e502c8eaa35096f9d0afbbaa3b97ca579a0a125f8b0a321` | Computation unchanged; docstring clarified |
| `ISCIENCE_REVISION_20260711/gee_validation/run_gee_primary_sensitivity_20260719_203058.py` | Same relative path | `7bfef1126ee62ba76e2b5f31fa0ed037a84bc189fe17dacd99ffce15cde8f491` | Byte-identical |
| `ISCIENCE_REVISION_20260711/gee_validation/run_gee_structured_null_refinement_20260719_205905.py` | Same relative path | `f6b8a073bd5ea7a6a55da164ccdb844ec98a844d45db46207ff2e5d43743a140` | Byte-identical |
| `ISCIENCE_REVISION_20260711/supplement/build_reader_facing_tables_s2_s3_20260718_235342.py` | Same relative path | `696f0c7185db23599c0a0fd9ecacae38d55bb4ef9391dc7e8cf903004f00e191` | Portable shebang only |
| `ISCIENCE_REVISION_20260711/supplement/build_supplemental_legends_v7_20260719_211018.py` | Same relative path | `08bbea86c7d8afbd9d6d337315293d27b6a7743f62d09a245881ca72187e5131` | Portable shebang only; retained V8 helper dependency |
| `ISCIENCE_REVISION_20260711/supplement/build_supplemental_legends_v8_20260719_214546.py` | Same relative path | `3cb5c5568b5ceb037099b7b1d65f9f1b75b8824b8cde3a2f5d85c600bf4f687a` | Portable shebang only |
| `ISCIENCE_REVISION_20260711/supplement/build_table_s2_v6_sha256_expansion_20260719_215020.py` | Same relative path | `9b5a038ce7a1f861111c1c23e9ec2861bf5efa9b7d120e0d018deb8433bb3998` | Portable shebang only |
| `ISCIENCE_REVISION_20260711/supplement/build_tables_s2_v3_s3_v6_20260719_210611.py` | Same relative path | `8eaa7aea1cf6b9dc0b424af6e21ff8859cf0ed8c75e0a72b3b6172a18cfbfcbc` | Portable shebang only |
| `ISCIENCE_REVISION_20260711/supplement/build_tables_s2_v4_s3_v7_20260719_213408.py` | Same relative path | `1fa3cc46b837b2a69b69af2e426f8420b7ff76412389e170c74e2c522a5bf9b0` | Portable shebang only |
| `ISCIENCE_REVISION_20260711/supplement/build_tables_s2_v5_s3_v8_20260719_214045.py` | Same relative path | `ad7a69a9c66f90ff08b537157ee1bae311e1511cefc444d6fc639a4b66f222cc` | Portable shebang only |

The timestamped `public_release_inventory_*.tsv` generated for v1.0.2 records the SHA-256 digest and byte size of every released file, including all retained v1.0.1 dependencies.
