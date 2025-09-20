# Taken the set of all proteomes of the analyzed Staphylococcus species, create a MMseqs2 database for further analysis.
mmseqs createdb \
     /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/new_extraction/only_proteomes/*.faa \
     staphylococci_proteome/all_proteomes_DB