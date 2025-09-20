
def read_mmseqs_homology(path): 
    with open(path, 'r') as f: 
        lines = []
        for line in f:
            split_line = line.split('\t')
            lines.append(split_line)
    return lines

def read_json(path):
    import json
    with open(path, 'r') as f:
        return json.load(f)
    
def adapt_mmseqs(lines, gff_name_to_locus, locus_to_organism):
    adapted_lines = []
    for line in lines:
        # get locus_tag
        locus_tags = gff_name_to_locus.get(line[1], None)
        for locus_tag in locus_tags:
            organism = locus_to_organism.get(locus_tag.split('_')[0], None)
            if locus_tag and organism:
                adapted_lines.append([line[0], line[1], locus_tag, organism, line[2], line[3], line[4], line[5], line[6], line[7], line[8]])
    return adapted_lines

def write_adapted_mmseqs(adapted_lines, path):
    with open(path, 'w') as f:
        for line in adapted_lines:
            f.write('\t'.join(line) + '\n')

def main():
    mmseqs_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_results_staphylococci_with_cntA_new.tsv"
    gff_name_to_locus_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/new_extraction/only_gffs/gff_name_to_locus.json"
    locus_to_organism_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/new_extraction/only_gffs/locus_to_organism.json"
    adapted_mmseqs_path = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/homology_analysis_MMSEQS/homology_results_staphyloccus_adapted_with_cntA_new.tsv"
    mmseqs = read_mmseqs_homology(mmseqs_path)
    gff_name_to_locus = read_json(gff_name_to_locus_path)
    locus_to_organism = read_json(locus_to_organism_path)
    adapted_mmseqs = adapt_mmseqs(mmseqs, gff_name_to_locus, locus_to_organism)
    write_adapted_mmseqs(adapted_mmseqs, adapted_mmseqs_path) 

if __name__ == '__main__':
    main()