for file in /supplementary/accession_codes/*_with_assembly_completeness.csv; do
  # get the first column of the file and remove the header
  cut -f1 "$file" | tail -n +2 > "${file%.csv}_genome_ids.txt"
  # for each line in the file, run the following command
  for line in $(cat "${file%.csv}_genome_ids.txt"); do
    # get the metadata for the genome and save it to a file
    datasets summary genome accession "$line" > "${file%.csv}_${line}_metadata.json"
  done
done