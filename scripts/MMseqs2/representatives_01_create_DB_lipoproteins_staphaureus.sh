# Given a set of protein sequences in FASTA format, create a MMseqs2 database for further analysis.
# This case, reference: Staphylococcus aureus lipoproteins of interest
mmseqs createdb \
     /ceph/ibmi/it/projects/CMFI/Data/saureus_proteins/fasta_proteins_saureus.fa \
     /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/lipoproteins_interest/lipoproteins_DB