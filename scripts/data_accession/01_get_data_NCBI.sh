file="supplementary/accession_codes/adapt_no_redundant.csv"
cut -f1 "$file" | tail -n +2 > "${file%.csv}_genome_ids.txt"
datasets download genome accession \
    --inputfile "${file%.csv}_genome_ids.txt" \
    --filename staphyloccoci_genomes.zip \
    --include genome,gbff,gff3,protein