from ast import parse
import pandas as pd
import Bio.SeqIO as SeqIO
def read_mmseqs(path, IDs_to_use):
    """
    Reads a mmseqs2 output file and returns a list of the protein IDs and the original sequence ID.

    Args:
        path (str): Path to the mmseqs2 output file.

    Returns:
        a list of tuples (protein_id, original_sequence_id)
    """
    map_to_hits = {}
    
    mmseqs_df = pd.read_csv(path, sep="\t", header=None)
    # set header of the dataframe
    mmseqs_df.columns = ["lipoprotein", "hit", "hit_locus_tag", "hit_genome_id", "bits", "identity", "evalue", "start_alignment", "end_alignment", "length_alignment", "differences"]
    print(mmseqs_df.head())
    # set identity to float
    mmseqs_df["identity"] = mmseqs_df["identity"].astype(float)
    # filter for identity > 50
    mmseqs_df = mmseqs_df[mmseqs_df["identity"] > 0.5]
    # filter for hit_genome_id in IDs_to_use
    mmseqs_df = mmseqs_df[mmseqs_df["hit_genome_id"].isin(IDs_to_use)]
    print(f"Number of hits with identity > 50: {len(mmseqs_df)}")
    # sort by identity descending
    mmseqs_df = mmseqs_df.sort_values(by=["identity"], ascending=False)
    # drop duplicates based on the hit_genome_id column, keeping the first occurrence
    mmseqs_df = mmseqs_df.drop_duplicates(subset=["lipoprotein", "hit_genome_id"], keep="first")
    # keep only htsA and sirA hits
    mmseqs_df = mmseqs_df[mmseqs_df["lipoprotein"].str.contains("HtsA|SirA|CntA")]
    print(f"Number of HtsA and SirA hits: {len(mmseqs_df)}")
    print(mmseqs_df.head())
    # iterate over the dataframe and append the protein_id
    for index, row in mmseqs_df.iterrows():
        if row["lipoprotein"] in map_to_hits:
            map_to_hits[row["lipoprotein"]].append((row["hit"], row["hit_locus_tag"], row["hit_genome_id"]))
        else:
            map_to_hits[row["lipoprotein"]] = [(row["hit"], row["hit_locus_tag"], row["hit_genome_id"])]
    return map_to_hits

def write_fasta(sequences, output_path):
    pass

def parse_proteome(path):
    proteome_dict = {}
    for record in SeqIO.parse(path, "fasta"):
        proteome_dict[record.id] = str(record.seq)
    return proteome_dict

def parse_metadata(path):
    df_metadata = pd.read_csv(path, sep="\t", header=0)
    df_metadata_sub = df_metadata[["accession", "gtdb_taxonomy"]]
    df_metadata_sub["species"] = df_metadata_sub["gtdb_taxonomy"].apply(lambda x: x.split("s__")[-1])
    return df_metadata_sub.set_index("accession")["species"].to_dict()

def parse_ids(path):
    df_ids = pd.read_csv(path, header=0)
    return df_ids["accession_code"].tolist()
def main():
    mmseqs_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_results_staphyloccus_adapted_with_cntA_new.tsv"
    output_fasta_path = "path/to/output.fasta"
    proteome_dir = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/new_extraction/only_proteomes/"
    path_to_ids = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/repo/StaphyDiversity/figures/Figure2/input/contigs_number.csv"
    ids_to_use = parse_ids(path_to_ids)
    sequences = read_mmseqs(mmseqs_path, ids_to_use)
    genome_to_species = parse_metadata("/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/gtdb-search_no_redundant.tsv")
    # get length of HtsA and SirA arrays
    cnta_length = len(sequences.get("CntA", []))
    htsa_length = len(sequences.get("HtsA", []))
    sira_length = len(sequences.get("SirA", []))
    print(f"Number of HtsA hits: {htsa_length}")
    print(f"Number of SirA hits: {sira_length}")
    print(f"Number of CntA hits: {cnta_length}")
    # Agress with analysis of heatmap. 

    # For each key of the sequences dictionary, 
    # retrieve the corresponding sequences from the proteome files in proteome_dir
    # and write them to the output_fasta_path in fasta format.
    for lipoprotein, hits in sequences.items():
        with open(f"{lipoprotein}_hits.fasta", "w") as f:
            for hit, locus_tag, genome_id in hits:
                proteome_path = f"{proteome_dir}/{genome_id}.faa"
                parsed_proteome = parse_proteome(proteome_path)
                species_name = genome_to_species[genome_id]
                # replace spaces in species_name with underscores
                species_name = species_name.replace(" ", "_")
                # get the sequence for the locus_tag from the parsed_proteome
                sequence = parsed_proteome.get(hit, None)
                if sequence:
                    # adapt sequence to break after 80 characters
                    sequence = '\n'.join(sequence[i:i+80] for i in range(0, len(sequence), 80))
                    f.write(f">{locus_tag}|{genome_id}|{species_name}\n{sequence}\n")
                else:
                    print(f"Warning: {hit} not found in {proteome_path}")
               


if __name__ == "__main__":
    main()