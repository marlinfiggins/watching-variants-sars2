# Useful scripts for analyzing variant data

## Prepping sequence counts from curated metadata

We can reduce the pruned metadata file generated in `../data` using the `prep_sequence_counts` script.
This takes in a configuration file specifying the countries and divisions to analyze, the date range to analyze, and the mapping for clades and lineages to analyze.
An example config can be found at `insert-config-path-here.yaml`

```shell
python prep_sequence_counts.py --config <path to config> --metadata-path <path to metadata> --export-path <path to export .tsv>
```

This will produce a .tsv file containing columns `date`, `location`, `variant`, `sequences` where each row corresponds to the number of sequences collected of the given variant on that day in the given location.
