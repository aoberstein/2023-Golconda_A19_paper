#!/bin/bash

# files=$(find . -iname '*metadata_GSE*.tab')
# files=(`find . -iname '*metadata_GSE*.tab'`)
files=(`find . -iname 'metadata_phred_strand_rep_collapsed.tab'`)

# echo ${files[@]}

## grab header from first metadata file in array
rm metadata_merged.tab
for f in ${files[0]}
    do head -n 1 $f > metadata_merged.tab
done

## write all metadata files, skipping the header, to "metadata_merged.tab"
for f in ${files[@]}
    do tail -n +2 $f >> metadata_merged.tab
done
