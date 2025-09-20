
# Set analysis folder
folder="/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/antiSMASH"
data="/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/ncbi_dataset/data/"
# Create results folder if it doesn't exist
mkdir -p $folder/results/
for file in $data/**/*.gbff; do
    echo $file
    # Get name of parent folder
    # This is the same as the accession code of the genome
    parent=$(basename $(dirname $file))
    # Check if parent folder exists in results, if yes, skip
    # Helps if the script crashes and you need to restart
    if [ -d "$folder/results/$parent" ]; then
        echo "Folder exists, skipping"
        continue
    fi
    # antismash v7.1 was used
    antismash $file \
    --output-dir $folder/results/$parent \
    --output-basename "$parent" \
    --html-start-compact \
    --cc-mibig \
    --cb-general \
    -c 12 \
    --cb-subclusters \
    --cb-knownclusters \
    --fullhmmer \
    --tigrfam \
    --tfbs \
    --genefinding-tool none
done

## The same script was used to analyse the strains of the selected species