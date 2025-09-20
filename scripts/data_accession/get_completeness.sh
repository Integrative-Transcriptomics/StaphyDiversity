
for file in ../*.txt; do
# get filename without path
    filename=$(basename -- "$file")
    echo $filename
    # modify filename to include _with_assembly_completeness before .txt
    filename="${filename%.txt}_with_assembly_completeness.txt"
    
    # for each line of given file, run the following command
    while read line; do
        datasets summary genome accession $line --as-json-lines
    done < $file > $filename 
done