# Hubmap registration scripts

This repository contains utility functions, scripts, and notebooks for data and metadata registration of [laminlabs/hubmap](https://lamin.ai/laminlabs/hubmap).

The order that was used for the single-cell data is as follows:

1. Generate a metadata table using the API (`generate-single-cell-metadata-table.ipynb`).
   This results in a first metadata csv table.
2. Curate and register all metadata (`register-metadata.ipynb`).
   The results in a second curated metadata table.
3. Register the Artifacts using the URLs that are stored in the metadata table and associated the metadata records with them (`register-single-cell-rna.ipynb`).
4. Optionally register Vitessce configs.
