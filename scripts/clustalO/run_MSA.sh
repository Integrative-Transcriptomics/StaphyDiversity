clustalo -i /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/2025_09_27_alignments_lipoproteins/HtsA_hits.fasta\
    -o output_htsa.aln \
    --outfmt=clu \
    --force

clustalo -i /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/2025_09_27_alignments_lipoproteins/SirA_hits.fasta\
    -o output_sira.aln \
    --outfmt=clu \
    --force
clustalo -i /ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/2025_09_27_alignments_lipoproteins/CntA_hits.fasta\
    -o output_cnta.aln \
    --outfmt=clu \
    --force