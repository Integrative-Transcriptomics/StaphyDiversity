# Identification of the threshold for homologous sequences within Staphylococcus species
# Using the previously created database of Staphylococcus aureus lipoproteins of interest, we
# perform a search against itself to identify homologous sequences within the Staphylococcus genus.
# Based on non self-hits of the lipoproteins of sequences we identify the minimal threshold. 
mmseqs search \
     -a --alignment-mode 3  --max-seqs 1000000 \
          /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/lipoproteins_interest/lipoproteins_DB\
          /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/lipoproteins_interest/lipoproteins_DB\
    /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_results/homology_within_staphs \
    tmp \
    --threads 24 
  