import pandas as pd
import os
import argparse

def parse_gmt(gmt_path: str, id_dict_path: str, output_path: str = None) -> pd.DataFrame:
    """Parse a GMT file and map gene IDs to Ensembl IDs, optionally saving results."""
    
    # Load the ID conversion dictionary, ensuring exactly two columns
    gene_dict = pd.read_csv(id_dict_path, sep='\t', dtype=str)
    if gene_dict.shape[1] != 2:
        raise ValueError(f"Expected 2 columns, but found {gene_dict.shape[1]}.")

    # Rename columns dynamically by checking for 'ens' (case-insensitive) in one of the columns
    gene_dict = gene_dict.rename(columns=lambda col: 'Ensembl_ID' if 'ensembl' in col.lower() else 'Other_ID')

    # Clean up: remove NaNs, duplicates, strip whitespace, and set index to 'Other_ID'
    gene_dict = (gene_dict.dropna(subset=['Other_ID'])
                         .drop_duplicates(subset=['Other_ID'])
                         .assign(Other_ID=lambda df: df['Other_ID'].str.strip())
                         .set_index('Other_ID'))

    # Initialize a list to store pathway data and a set to track missing IDs
    pathway_data = []
    missing_ids = set()

    # Process the GMT file: extract pathway names, descriptions, and convert IDs
    with open(gmt_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # Skip lines with insufficient data
            
            pathway_name, description, *other_ids = parts
            other_ids = map(str.strip, other_ids)  # Clean up gene names

            for other_id in other_ids:
                if other_id in gene_dict.index:
                    ensembl_id = gene_dict.at[other_id, 'Ensembl_ID']
                    pathway_data.append([pathway_name, description, ensembl_id])
                else:
                    missing_ids.add(other_id)

    # Report number of missing IDs
    print(f"Total pathway genes not found in dictionary: {len(missing_ids)}")

    # Create a DataFrame from the collected pathway data
    pathway_df = pd.DataFrame(pathway_data, columns=['Pathway_Name', 'Description', 'Ensembl_ID'])

    # If no output path is provided, return the DataFrame without saving
    if output_path is None:
        return pathway_df

    # Create the output filename based on the GMT input file
    gmt_filename = os.path.basename(gmt_path).replace('.gmt', '')
    output_file = os.path.join(output_path, f'{gmt_filename}_parsed_ensembl.tsv')

    # Save the DataFrame as a TSV file
    pathway_df.to_csv(output_file, sep='\t', index=False)

    print(f"Results saved to {output_file}")

    return pathway_df

def parse_gmt_main():
    """Parse command-line arguments and execute the parse_gmt function."""
    
    parser = argparse.ArgumentParser(description="Parse GMT file and map gene IDs to Ensembl IDs.")
    parser.add_argument('-g', '--gmt_path', required=True, help="Path to the GMT file.")
    parser.add_argument('-i', '--id_dict_path', required=True, help="Path to the ID dictionary (tab-separated file).")
    parser.add_argument('-o', '--output_path', default=None, help="Directory to save the parsed TSV file.")
    
    args = parser.parse_args()

    # Execute the parse_gmt function with the parsed arguments
    parse_gmt(
        gmt_path=args.gmt_path,
        id_dict_path=args.id_dict_path,
        output_path=args.output_path
    )

if __name__ == "__main__":
    # parse_gmt_main()
    # Uncomment the previous line and comment the following lines to run from the command line

    # Direct invocation for IDE usage
    gmt_path = '../../../data/pathways/KEGG/KEGG_2022.gmt'
    id_dict_path = '../../../data/id_conversion/ENTREZID_to_ENSEMBL.txt'
    output_path = '../../../data/pathways/KEGG/'
    
    pathway_df = parse_gmt(
        gmt_path=gmt_path,
        id_dict_path=id_dict_path,
        output_path=output_path
    )
