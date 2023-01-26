# Data preparation

## Preparing metadata for sequence counts

To begin data preparation, we first download the Nextstrain-curated metadata TSV of the GISAID database.
This metadata file is used to construct the datasets in the sub-folders.

1. Nextstrain-curated metadata TSV from the GISAID database was downloaded using [nextstrain-cli](https://docs.nextstrain.org/projects/cli/en/stable/). Uncompressing and renaming this file resulted in `gisaid_metadata.tsv` via:
```
nextstrain remote download s3://nextstrain-ncov-private/metadata.tsv.gz
gzip -d metadata.tsv.gz -c > gisaid_metadata.tsv
```

2. The metadata file was pruned to only relevant columns using the `tsv-select` command from [tsv-utils](https://opensource.ebay.com/tsv-utils/):
```
tsv-select -H -f strain,date,country,division,QC_overall_status,Nextstrain_clade,Nextclade_pango gisaid_metadata.tsv > gisaid_metadata_pruned.tsv
```

These operations can be found in the script `./download_prune_metadata.sh`.

## Using open data

If you would like to reproduce this work and do not have access to the Nextstrain-curated metadata. Substituing `s3://nextstrain-ncov-private/metadata.tsv.gz` for `s3://nextstrain-data/files/ncov/open/metadata.tsv.gz` will instead use Nextstrain-curated open data from Genbank instead of GISAID data. Details on the various endpoints can be found on [Nextstrain remote inputs](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html).
