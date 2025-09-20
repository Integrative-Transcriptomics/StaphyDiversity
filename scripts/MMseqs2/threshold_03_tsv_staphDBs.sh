# Create a TSV file based on the results of the previous search against itself to identify homologous sequences within Staphylococcus species.
mmseqs createtsv \
     /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/lipoproteins_interest/lipoproteins_DB \
     /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/lipoproteins_interest/lipoproteins_DB \
    /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_results/homology_within_staphs \
    /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_within_staphs.tsv
