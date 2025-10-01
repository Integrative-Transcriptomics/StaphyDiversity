# A landscape of metallophore synthesis and uptake potential of the genus Staphylococcus
## Scripts 

The following scripts are classified with respect to the tools used.

- antiSMASH (v7.1.0): It was used locally for the identification and annotation of BGCs. A further script was used to parse all produced JSON files and summarize the results into one file. The script only accounts for BGCs labelled as metallophores or siderophores. 
- bigscape (v2): The tool BiG-SCAPE was used to query for one specific BGC in multiple results. This was only used to identify the presence of a known BGC in the results of S. lugdunensis at the strain level. 
- chewBBACA (v3.3.10): The MLST tool chewBBACA was used to produce a distance matrix of the samples analyzed. The desktop tool MEGA (v11.0.11) was used afterwards to reconstruct a phylogenetic tree using a NJ algorithm.
- clustalO (v1.2.4): For the alignemnt of multiple homologs, we used a script to extract the homologs from our MMseqs2 results. The FASTA file was then provided to clustal Omega, which aligned the genes. Based on the alignment and knowledge from the literature, we identified the conservation of the binding sites based on the sequence originating from S. aureus.
- data accession: Scripts used for the data collection. 
- MMseqs2 (v15, commit ad6dfc): Identification of homologs based on MMseqs2. It shows the scripts used for the identification of homologs across the 77 representatives (representatives_), the four strain analyses (strains_) but also for the computation of the homology threshold (threshold_)