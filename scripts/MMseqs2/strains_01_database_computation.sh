# Create DB based on proteomes
mmseqs createdb /beegfs/work/workspace/ws/Analaysis_Staph_SIDS-0/Hominis/data/ncbi_dataset/data/**/protein*.fasta \
    /beegfs/work/workspace/ws/Analaysis_Staph_SIDS-0/Hominis/analysis/mmseqs/dbmmseqs/HominisDB
# Index computation 
mmseqs createindex /beegfs/work/workspace/ws/Analaysis_Staph_SIDS-0/Hominis/analysis/mmseqs/dbmmseqs/HominisDB tmp