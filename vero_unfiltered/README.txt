R2 reads were discarded. Reads from Lane1 and Lane2 (L001, L002) were concatenated and saved into fastq/.

The sgRNA-library was taken from the 'Vervet_KO_lib_with_new_SeqID.xlsx' file. As can be seen in the results, I used gene symbols where available and fell back to ENSEMBL ids, where symbols were not available. It should be easy to merge back information like targeted-domains to the results.

Different to the other Vero analysis, here I did not filter any sgRNAs due to multiple targets. Yet I merged targets if they were targeted by the same sgRNA(s), giving them names starting with multiKO_
