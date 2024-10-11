import requests
import pandas as pd

# Fetch the data from the URL
url = "https://rest.kegg.jp/conv/ncbi-geneid/hsa"
response = requests.get(url)

# Split the response text into lines and columns
data = [line.split('\t') for line in response.text.strip().split('\n')]

# Create a DataFrame
df = pd.DataFrame(data, columns=["KEGG_ID", "Entrez_ID"])

# Remove the 'ncbi-geneid:' prefix from the Entrez_ID column
df["Entrez_ID"] = df["Entrez_ID"].str.replace("ncbi-geneid:", "", regex=False)

# Specify the output path
output_path = '../../data/IDs/KEGG_ID_Entrez.txt'

# Save the DataFrame to a tab-separated text file
df.to_csv(output_path, sep='\t', index=False)

# Optionally, print the DataFrame
print(df)
