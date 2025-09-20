import json
import os
import pandas as pd
def read_json(path):
    print(path)
    with open(path, "r") as f:
        return json.load(f)
    
genes_BGCs = {
    # "BGC0002487": ['SAV2462', 'SAV2463', 'SAV2464', 'SAV2465', 'SAV2466', 'SAV2467', 'SAV2468', 'SAV2469', 'SAV2470'],
    # Staphylopine genes not involved in uptake 
    # cntE, cntM, cntL, cntK 
    "BGC0002487": ['SAV2462', 'SAV2468', 'SAV2469', 'SAV2470'],
    # Stapphyloferrin A: sfaC, sfaB, sfaA, sfaD
    "BGC0000944": ['SAOUHSC_02433', 'SAOUHSC_02434', 'SAOUHSC_02435', 'SAOUHSC_02436'],
    # Staphyloferrin B: sbnA, sbnB, sbnC, sbnD, sbnE, sbnF, sbnG, sbnH, sbnI
    "BGC0000943": ['NCTC8325_00065', 'NCTC8325_00066', 'NCTC8325_00067', 'NCTC8325_00068', 'NCTC8325_00069', 'NCTC8325_00070', 'NCTC8325_00071', 'NCTC8325_00072', 'NCTC8325_00073']
}


def summarize_json(json_file):
    all_records = json_file["records"]
    interesting_records = {}
    for index_record, record in enumerate(all_records): 
        save_records = {}
        for i,area in enumerate(record["areas"], 1):
            save_record = {}
            # if any product contains "ophore" in the name, save the record
            for product in area["products"]:
                if "ophore" in product:
                    save_record["record_key"] = f"record_{index_record}_{i}"
                    save_record["index"] = i
                    save_record["product"] = product
                    save_records[i] = save_record
        interesting_records[index_record] = save_records
        gathered_interesting_indices = [record_unit["index"] for record_unit in save_records.values()]
        if not "antismash.modules.cluster_compare" in record["modules"]:
            # No results based on cluster compare
            continue
        hits_clusterblast = record["modules"]["antismash.modules.cluster_compare"]
        hits_cluster_compare_per_region = hits_clusterblast["db_results"]["MIBiG"]["by_region"]
       

        for region in hits_cluster_compare_per_region.keys():
            hits_with_score = hits_cluster_compare_per_region[region]["RegionToRegion_RiQ"]["scores_by_region"]
            # From {key:score} to sorted list of tuples
            hits_with_score = sorted(hits_with_score.items(), key=lambda x: x[1], reverse=True)[:5]
            filtered_hits_with_score = [hit[0] for hit in hits_with_score if ("BGC0002487" in hit[0] or "BGC0000944" in hit[0] or "BGC0000943" in hit[0]) and (hit[1] > 0.3)]
            if not int(region) in gathered_interesting_indices:
                continue
            hits_region = hits_cluster_compare_per_region[region]["RegionToRegion_RiQ"]["hits"]
            if len(filtered_hits_with_score) == 0:
                continue
            else:
                region_number = int(region)
                found_hit = filtered_hits_with_score[0]
               
                BGC_ID = found_hit.split(":")[0]
                
                hit_data = hits_region[found_hit]
                genes_BGC = genes_BGCs[BGC_ID] if BGC_ID in genes_BGCs else []
                filtered_genes = [gene for gene in genes_BGC if gene in hit_data.keys() and hit_data[gene]["percent_identity"] > 0.3]
                missing_genes = [gene for gene in genes_BGC if gene not in filtered_genes]
                completeness = len(filtered_genes) / len(genes_BGC) if len(genes_BGC) > 0 else 0
                interesting_records[index_record][region_number]["clustercompare"] = {}
                interesting_records[index_record][region_number]["clustercompare"]["BGC_ID"] = BGC_ID
                interesting_records[index_record][region_number]["clustercompare"]["genes"] = filtered_genes
                interesting_records[index_record][region_number]["clustercompare"]["missing_genes"] = missing_genes
                interesting_records[index_record][region_number]["clustercompare"]["completeness"] = completeness*100
                
               

    flattened_records = []
    for record_key, record_value in interesting_records.items():
        for record_unit in record_value.values():
            if not "clustercompare" in record_unit:
                record_unit["clustercompare"] = "No hits"
            flattened_records.append(record_unit)
    return flattened_records

def analyze_directory(path, metadata_path):
    parsed_metadata = pd.read_csv(metadata_path, sep="\t")
    dict_metadata = parsed_metadata[["accession", "ncbi_organism_name"]].set_index("accession").to_dict()["ncbi_organism_name"]
    all_samples = {}
    for sample in os.listdir(path):
        if "json" in sample:
            continue
        if not os.path.exists(f"{path}/{sample}/{sample}.json"):
            continue
        full_data = summarize_json(read_json(f"{path}/{sample}/{sample}.json"))

        unknown_hits = [record for record in full_data if record["clustercompare"] == "No hits"]
        known_hits = [record for record in full_data if record["clustercompare"] != "No hits"]
        Sfa = [record["clustercompare"]["completeness"] for record in known_hits if "BGC0000944" == record["clustercompare"]["BGC_ID"]]
        Sfb = [record["clustercompare"]["completeness"] for record in known_hits if "BGC0000943" == record["clustercompare"]["BGC_ID"]]
        Stp = [record["clustercompare"]["completeness"] for record in known_hits if "BGC0002487" == record["clustercompare"]["BGC_ID"]]
        all_samples[sample] = {
            "genome": dict_metadata[sample],
            "data": full_data, 
            "number_of_records": len(full_data),
            "unknown" : len(unknown_hits),
            "staphyloferrin A" : Sfa if len(Sfa) > 1 else Sfa[0] if len(Sfa) == 1 else 0,
            "staphyloferrin B" : Sfb if len(Sfb) > 1 else Sfb[0] if len(Sfb) == 1 else 0,
            "staphylopine" : Stp if len(Stp) > 1 else Stp[0] if len(Stp) == 1 else 0
        }

    return all_samples

def adapt_output(path, output_path, metadata): 
    all_samples = analyze_directory(path, metadata)
    with open(f"{output_path}/summary_with_clusterblast_with_percentage_clustercompare.json", "w") as f:
        json.dump(all_samples, f)
    # create a dataframe from the json
    df = pd.DataFrame.from_dict(all_samples, orient="index")
    # remove column data
    df = df.drop(columns=["data"])
    df.to_csv(f"{output_path}/summary_with_clusterblast_with_percentage_clustercompare.csv", sep=";")

def main():
    capitis_input = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/capitis/results/results"
    capitis_output = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/final_results/antismash/capitis"
    capitis_metadata = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/data/S_capitis_all.tsv"
    adapt_output(capitis_input, capitis_output, capitis_metadata)

    hominis_input = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/hominis/results"
    hominis_output = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/final_results/antismash/hominis"
    hominis_metadata = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/data/S_hominis_all.tsv"
    adapt_output(hominis_input, hominis_output, hominis_metadata)

    lugdunensis_input = "/ceph/ibmi/it/projects/CMFI/Strains_Staphyloccocus_Analysis/analysis/more_analysis/lugdunensis/test_antismash/antismash/results"
    lugdunensis_output = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/final_results/antismash/lugdunensis"
    lugdunensis_metadata = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/data/S_lugdunensis_all.tsv"
    adapt_output(lugdunensis_input, lugdunensis_output, lugdunensis_metadata)

    epi_input = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/epidermidis/results"
    epi_output = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/final_results/antismash/epidermidis"
    epi_metadata = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/04_11_2024_strain_analysis/data/S_epidermidis_all.tsv"
    adapt_output(epi_input, epi_output, epi_metadata)

    all_input = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/analysis/antiSMASH/results"
    all_output = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/results/final_results/antismash/all_species"
    all_metadata = "/ceph/ibmi/it/projects/CMFI/Staphyloccus_Analysis/data/gtdb-search.tsv"
    adapt_output(all_input, all_output, all_metadata)



    

main()