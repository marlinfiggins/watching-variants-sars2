# Data preparation

## Creating data set focused on rise of XBB.1.5. in the US.

This dataset can be created with the command

```shell
python3 ../../scripts/prep_sequence_counts.py  --config config.yaml \ 
--metadata-path ../../data/gisaid_metadata_pruned.tsv \
--export-path ./us-xbb15-seq_counts.tsv
```

