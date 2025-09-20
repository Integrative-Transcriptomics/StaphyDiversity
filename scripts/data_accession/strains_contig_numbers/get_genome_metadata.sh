for file in /ceph/ibmi/it/projects/CMFI/Strains_Staphyloccocus_Analysis/data/S_*_all.tsv; do
  # get the first column of the file and remove the header
  cut -f1 "$file" | tail -n +2 > "${file%.tsv}_genome_ids.txt"
  # for each line in the file, run the following command
  for line in $(cat "${file%.tsv}_genome_ids.txt"); do
    # get the metadata for the genome and save it to a file
    datasets summary genome accession "$line" > "${file%.tsv}_${line}_metadata.json"
  done
done