# Function to parse a GMT file, including descriptions
def parse_gmt(gmt_file):
    gmt_data = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            pathway_id = parts[0]
            description = parts[1]
            genes = parts[2:]  # Genes are listed from the third element onward
            gmt_data[pathway_id] = {"description": description, "genes": genes}
    return gmt_data

# List of pathway identifiers to filter by
infection_pathways = [
    'hsa05166',  # Human T-cell leukemia virus 1 infection
    'hsa05170',  # Human immunodeficiency virus 1 infection
    'hsa05161',  # Hepatitis B
    'hsa05160',  # Hepatitis C
    'hsa05171',  # Coronavirus disease - COVID-19
    'hsa05164',  # Influenza A
    'hsa05162',  # Measles
    'hsa05168',  # Herpes simplex virus 1 infection
    'hsa05163',  # Human cytomegalovirus infection
    'hsa05167',  # Kaposi sarcoma-associated herpesvirus infection
    'hsa05169',  # Epstein-Barr virus infection
    'hsa05165',  # Human papillomavirus infection
    'hsa05110',  # Vibrio cholerae infection
    'hsa05120',  # Epithelial cell signaling in Helicobacter pylori infection
    'hsa05130',  # Pathogenic Escherichia coli infection
    'hsa05132',  # Salmonella infection
    'hsa05131',  # Shigellosis
    'hsa05135',  # Yersinia infection
    'hsa05133',  # Pertussis
    'hsa05134',  # Legionellosis
    'hsa05150',  # Staphylococcus aureus infection
    'hsa05152',  # Tuberculosis
    'hsa05100',  # Bacterial invasion of epithelial cells
    'hsa05146',  # Amoebiasis
    'hsa05144',  # Malaria
    'hsa05145',  # Toxoplasmosis
    'hsa05140',  # Leishmaniasis
    'hsa05142',  # Chagas disease
    'hsa05143'   # African trypanosomiasis
]

# Load and parse the GMT file
gmt_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez.gmt'
gmt_data = parse_gmt(gmt_file)

# Filter the GMT data by pathway identifiers
filtered_gmt = {pathway: data for pathway, data in gmt_data.items() if pathway in infection_pathways}

# Save the filtered pathways back to a new GMT file, including description
output_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_infection.gmt'
with open(output_file, 'w') as f:
    for pathway, data in filtered_gmt.items():
        description = data['description']
        genes = data['genes']
        f.write(f"{pathway}\t{description}\t" + "\t".join(genes) + "\n")

print(f"Filtered pathways saved to {output_file}")

