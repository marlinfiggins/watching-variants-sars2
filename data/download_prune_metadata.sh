nextstrain remote download s3://nextstrain-ncov-private/metadata.tsv.gz
gzip -d metadata.tsv.gz -c > gisaid_metadata.tsv
tsv-select -H -f strain,date,country,division,QC_overall_status,Nextstrain_clade,Nextclade_pango gisaid_metadata.tsv > gisaid_metadata_pruned.tsv
