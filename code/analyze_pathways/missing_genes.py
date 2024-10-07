import pandas as pd

def load_id_dict(ensembl_dict_path):
    gene_dict = pd.read_csv(ensembl_dict_path, sep='\t', header=0,
                             names=['Gene_Name', 'Ensembl_ID'],
                             dtype={'Gene_Name': str}
                            ).dropna().drop_duplicates(subset=['Gene_Name']).set_index('Gene_Name')
    return gene_dict

def parse_gmt_file(gmt_path, gene_dict):
    pathway_rows = []
    error_summary = []
    all_invalid_genes = set()

    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathway_name, description, *gene_names = parts
                gene_names = list(set(gene_names))  # Remove duplicates
                
                valid_genes = [name for name in gene_names if name in gene_dict.index]
                invalid_genes = set(gene_names) - set(valid_genes)
                
                if valid_genes:
                    ensembl_ids = gene_dict.loc[valid_genes, 'Ensembl_ID'].unique()
                    pathway_rows.extend([[pathway_name, description, ensembl_id] for ensembl_id in ensembl_ids])

                error_summary.append({
                    'Pathway_Name': pathway_name,
                    'Description': description,
                    'Total_Genes': len(gene_names),
                    'Not_Converted_Count': len(invalid_genes),
                    'Not_Converted_Genes': list(invalid_genes)
                })
                
                all_invalid_genes.update(invalid_genes)

    valid_df = pd.DataFrame(pathway_rows, columns=['Pathway_Name', 'Description', 'Ensembl_ID'])
    error_df = pd.DataFrame(error_summary)
    invalid_gene_list = list(all_invalid_genes)

    return valid_df, error_df, invalid_gene_list

# Paths to input files
gmt_path = "../../data/pathways/KEGG/KEGG_2022_entrez.gmt"
ensembl_dict_path = "../../data/IDs/entrez2ensembl_kegg.txt"

# Load the gene dictionary and parse the GMT file
gene_dict = load_id_dict(ensembl_dict_path)
valid_df, error_df, invalid_gene_list = parse_gmt_file(gmt_path, gene_dict)

# Save valid_df to a CSV file
output_path = "../../data/pathways/KEGG/KEGG_2022_ensembl.csv"
valid_df.to_csv(output_path, index=False)

# Optionally, save the error_df and invalid gene list as well
error_output_path =  "../../data/pathways/KEGG/KEGG_2022_errors.csv"
error_df.to_csv(error_output_path, index=False)


print("Data saved successfully!")

    
    