import pandas as pd
def read_json(file_path):
    import json
    with open(file_path, 'r') as file:
        json_data = json.load(file)
        if "reports" not in json_data:
            # extract id from the file name
            # get the file name without the directory and extension
            import os
            file_name = os.path.basename(file_path)
            file_name_without_ext = os.path.splitext(file_name)[0]
            # get id between _all_ and _metadata
            id_start = file_name_without_ext.find('_all_') + len('_all_')
            id_end = file_name_without_ext.find('_metadata')
            genome_id = file_name_without_ext[id_start:id_end]
            return genome_id
        return json_data.get("reports")[0]
def get_assembly_stats(json_file): 
    """
    Extracts the assembly statistics from a JSON file.
    
    Args:
        json_file (str): Path to the JSON file.
        
    Returns:
        dict: Assembly statistics.
    """
    data = read_json(json_file)
    if type(data) is str:
        return {"accession": data}
    accession = data.get('accession')

    assembly_stats = data.get('assembly_stats')
    combine_stats = assembly_stats
    combine_stats["accession"] = accession

    return combine_stats


def main():
    import os
    import glob

    # Path to the directory containing the JSON files
    # json_dir = '/ceph/ibmi/it/projects/CMFI/Strains_Staphyloccocus_Analysis/data/capitis/metadata'
    # json_dir = '/ceph/ibmi/it/projects/CMFI/Strains_Staphyloccocus_Analysis/data/hominis/metadata'
    # json_dir = '/ceph/ibmi/it/projects/CMFI/Strains_Staphyloccocus_Analysis/data/epidermidis/metadata'
    json_dir = '/ceph/ibmi/it/projects/CMFI/Strains_Staphyloccocus_Analysis/data/lugdunensis/metadata'

    # Get all JSON files in the directory
    json_files = glob.glob(os.path.join(json_dir, '*.json'))

    # Extract assembly statistics from each JSON file
    assembly_stats = [get_assembly_stats(json_file) for json_file in json_files]
    data_frame = pd.json_normalize(assembly_stats)

    print(data_frame)
    # Save the DataFrame to a CSV file
    output_csv = os.path.join(json_dir, 'assembly_stats_pandas.csv')
    data_frame.to_csv(output_csv, index=False)


if __name__ == "__main__":
    main()