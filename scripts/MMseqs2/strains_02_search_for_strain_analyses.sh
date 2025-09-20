# Same approach as for the representatives, but a slightly different command. 
# Adapted to be used on HPCCLuster
# Same script was used for all strain analysis, here Hominis as and example
mmseqs easy-search \
    -a --alignment-mode 3  --max-seqs 1000000 \
    --cov-mode 2 -c 0.35 \
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov"\
    /beegfs/work/workspace/ws/Analaysis_Staph_SIDS-0/general_data/lipoproteins_StaphAureus_with_cnta.fasta \
    /beegfs/work/workspace/ws/Analaysis_Staph_SIDS-0/Hominis/analysis/mmseqs/dbmmseqs/HominisDB \
    /beegfs/work/workspace/ws/Analaysis_Staph_SIDS-0/Hominis/analysis/mmseqs/alnResult_last.m8 tmp