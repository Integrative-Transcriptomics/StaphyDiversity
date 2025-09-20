# Based on both created databases, perform a search to identify homologous sequences.
# Here, we use alignment mode 3 (more sensitive) and allow a maximum of 1,000,000 sequences to be reported per query.
# The results will be stored in the specified output file with the prefix "homology_DB_new_search_alignmentmode3".
mmseqs search \
     -a --alignment-mode 3  --max-seqs 1000000 \
     /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/lipoproteins_interest/lipoproteins_DB \
     /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/staphylococci_proteome/all_proteomes_DB\
    /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_results/homology_DB_new_search_alignmentmode3 \
    tmp \
    --threads 24 
  