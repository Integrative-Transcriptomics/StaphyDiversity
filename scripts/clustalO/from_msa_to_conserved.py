from Bio import AlignIO

def parse_alignment(path):
    alignment = AlignIO.read(path, "clustal")
    return alignment

def reduce_alphabet(aa, reduced_alphabet):
    translate_dictionary = {}
    new_sequence = ""
    match reduced_alphabet:
        case "Murphy8":
            translate_dictionary = {
                "L": "L",
                "V": "L",
                "I": "L",
                "M": "L",
                "C": "L",
                "A": "A",
                "G": "A",
                "S": "S",
                "T": "S",
                "F": "F",
                "Y": "F",
                "W": "F",
                "K": "K",
                "R": "K", 
                "E": "E",
                "D": "E",
                "N": "E",
                "Q": "E",
            }
        case "Murphy10": 
            translate_dictionary = {
                "L": "L",
                "V": "L",
                "I": "L",
                "M": "L",
                "S": "S",
                "T": "S",
                "F": "F",
                "Y": "F",
                "W": "F",
                "K": "K",
                "R": "K", 
                "E": "E",
                "D": "E",
                "N": "E",
                "Q": "E",
            }
        case "PC5": 
            translate_dictionary = translate_dictionary = {
                "L": "A",
                "V": "A",
                "I": "A",
                "C": "T",
                "A": "T",
                "G": "T",
                "S": "T",
                "F": "R",
                "Y": "R",
                "W": "R",
                "H": "R",
                "K": "C",
                "R": "C", 
                "E": "C",
                "D": "C",
                "T": "D",
                "P": "D",
                "M": "D",
                "N": "D",
                "Q": "D",
            }
        case "Murphy4": 
            translate_dictionary = translate_dictionary = {
                "L": "L",
                "V": "L",
                "I": "L",
                "M": "L",
                "C": "L",
                "A": "A",
                "G": "A",
                "S": "A",
                "T": "A",
                "P": "A",
                "F": "F",
                "Y": "F",
                "W": "F",
                "H": "E",
                "K": "E",
                "R": "E", 
                "E": "E",
                "D": "E",
                "N": "E",
                "Q": "E",
            }
    mod_aa = str(aa)
    if mod_aa in translate_dictionary:
        new_aa = translate_dictionary[mod_aa]
    else:
        new_aa = mod_aa
    return new_aa

HTSA_BINDING_SITES = [
    {"position": 86, "amino_acid": "R"},
    {"position": 104, "amino_acid": "R"},
    {"position": 126, "amino_acid": "R"},
    {"position": 203, "amino_acid": "K"},
    {"position": 209, "amino_acid": "H"},
    {"position": 239, "amino_acid": "Y"},
    {"position": 299, "amino_acid": "R"},
    {"position": 304, "amino_acid": "R"},
    {"position": 306, "amino_acid": "R"}
]
SIRA_BINDING_SITES = [
    {"position": 81, "amino_acid": "W"},
    {"position": 125, "amino_acid": "R"},
    {"position": 144, "amino_acid": "T"},
    {"position": 201, "amino_acid": "R"},
    {"position": 206, "amino_acid": "R"},
    {"position": 208, "amino_acid": "Y"},
    {"position": 304, "amino_acid": "N"},

]
CNTA_BINDING_SITES = [
    {"position": 27+25, "amino_acid": "Y"},
    {"position": 103+25, "amino_acid": "W"},
    {"position": 140+25, "amino_acid": "R"},
    {"position": 225+25, "amino_acid": "R"},
    {"position": 393+25, "amino_acid": "R"},
    {"position": 406+25, "amino_acid": "W"},
    {"position": 410+25, "amino_acid": "Y"},
    {"position": 423+25, "amino_acid": "N"},
    {"position": 497+25, "amino_acid": "Y"},

]

def get_position_MSA(sequence, position, ref_aa):
    subsequence_until_position = sequence[:position]
    gaps = subsequence_until_position.count("-")
    current_position = position + gaps -1
    aa_at_position = sequence[current_position].upper()
    while aa_at_position != ref_aa:
        subsequence_until_position = sequence[:current_position]
        gaps = subsequence_until_position.count("-")
        current_position = position + gaps -1
        aa_at_position = sequence[current_position].upper()
        if aa_at_position == "-":
            current_position += 1
   
    return current_position


def evaluate_alignment(alignment_path, lipoprotein, ID_aureus):
    alignment = parse_alignment(alignment_path)
    print(f"Number of sequences in alignment: {len(alignment)}")
    print(f"Alignment length: {alignment.get_alignment_length()}")
    # Get the sequence for S. aureus
    sequence_aureus = None
    for record in alignment:
        if ID_aureus in record.id:
            sequence_aureus = record.seq
            break

    binding_sites_in_MSA = {

    }
    reference_bindings = HTSA_BINDING_SITES if lipoprotein == "HtsA" else SIRA_BINDING_SITES if lipoprotein == "SirA" else CNTA_BINDING_SITES
    # get position of each binding site in the MSA
    for binding_site in reference_bindings:
        position = binding_site["position"]
        amino_acid = binding_site["amino_acid"]
        msa_position = get_position_MSA(sequence_aureus, position, amino_acid)
        print(f"Binding site at position {position} ({amino_acid}) in reference is at position {msa_position} in MSA")
        binding_sites_in_MSA[msa_position] = {"position": position, "amino_acid": amino_acid}
    print("Binding sites in MSA positions:", binding_sites_in_MSA)
    evaluation_by_position = {}
    # evaluate each sequence in the alignment
    for record in alignment:
        # if ID_aureus in record.id:
        #     continue
        for msa_position, msa_data in binding_sites_in_MSA.items():
            query_aa = record.seq[msa_position]
            id_evaluation = f"{msa_data["position"]}-{msa_data["amino_acid"]}"
            if id_evaluation not in evaluation_by_position:
                evaluation_by_position[id_evaluation] = {"identical":0, "modified_identical":0}
            data_position = evaluation_by_position.get(id_evaluation)
            if query_aa == msa_data["amino_acid"]:
                identical_positions = data_position.get("identical", 0)
                identical_positions += 1
                evaluation_by_position[id_evaluation]["identical"] = identical_positions
            if reduce_alphabet(query_aa, "Murphy10") == reduce_alphabet(msa_data["amino_acid"], "Murphy10"):
                mod_identical_positions = data_position.get("modified_identical", 0)
                mod_identical_positions += 1
                evaluation_by_position[id_evaluation]["modified_identical"] = mod_identical_positions
    print("Evaluation by position:", evaluation_by_position)

    #create table with the results
    total_sequences = len(alignment)  
    with open(f"binding_sites_evaluation_{lipoprotein}.csv", "w") as f:
        f.write("Position,Amino Acid,Identical,Modified Identical,Total Sequences,Identical %,Modified Identical %\n")
        for msa_position, msa_data in binding_sites_in_MSA.items():
            id_evaluation = f"{msa_data['position']}-{msa_data['amino_acid']}"
            data_position = evaluation_by_position.get(id_evaluation, {"identical":0, "modified_identical":0})
            identical = data_position.get("identical", 0)
            modified_identical = data_position.get("modified_identical", 0)
            identical_perc = (identical / total_sequences) * 100
            modified_identical_perc = (modified_identical / total_sequences) * 100
            f.write(f"{msa_data['position']},{msa_data['amino_acid']},{identical},{modified_identical},{total_sequences},{identical_perc:.2f},{modified_identical_perc:.2f}\n")

def main(): 
    ID_aureus = "GCF_001027105.1"
    alignment_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/2025_09_27_alignments_lipoproteins/output_htsa.aln"
    evaluate_alignment(alignment_path, "HtsA", ID_aureus)
    alignment_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/2025_09_27_alignments_lipoproteins/output_sira.aln"
    evaluate_alignment(alignment_path, "SirA", ID_aureus) 
    alignment_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/2025_09_27_alignments_lipoproteins/output_cnta.aln"
    evaluate_alignment(alignment_path, "CntA", ID_aureus)


if __name__ == "__main__":
    main()