import pandas as pd
import os
import re
import argparse

def parse_gff(gff_file, gene_name):
    """
    Parse GFF file to count occurrences of a specific gene where the gbkey is 'Gene'.
    
    Parameters:
        gff_file (str): Path to the GFF file.
        gene_name (str): Name of the gene to count.
        
    Returns:
        int: Count of gene occurrences where the gbkey is 'Gene'.
    """
    try:
        df = pd.read_csv(gff_file, sep="\t", comment='#', header=None,
                         names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
        # Match only rows where the attribute contains 'gbkey=Gene' and handle NA values
        df_genes = df[df['attribute'].str.contains('gbkey=Gene', na=False)]
        # Compile a regex to match the gene name in the attribute column with 'gene='
        gene_regex = re.compile(r'\bgene=' + re.escape(gene_name) + r'\b', re.IGNORECASE)
        gene_count = df_genes[df_genes['attribute'].apply(lambda x: bool(gene_regex.search(x)) if pd.notna(x) else False)].shape[0]
        return gene_count
    except FileNotFoundError:
        print(f"Error: File '{gff_file}' not found.")
        return 0
    except pd.errors.EmptyDataError:
        print(f"Error: File '{gff_file}' is empty or not properly formatted.")
        return 0

def main(output_path):
    genes = ["PRM1", "PRM2", "TNP1", "TNP2"]
    species = ["Homo_sapiens", "Mus_musculus", "Bos_taurus", "Sorex_araneus", "Zalophus_californianus", 
               "Tursiops_truncatus", "Physeter_catodon", "Balaenoptera_acutorostrata", 
               "Balaenoptera_musculus", "Loxodonta_africana", "Ornithorhynchus_anatinus"]
    data = {}

    for gene in genes:
        counts = []
        for sp in species:
            gff_filename = f"{sp}.gff"
            count = parse_gff(gff_filename, gene)
            counts.append(count)
        data[gene] = counts

    df_counts = pd.DataFrame(data, index=species)
    df_counts.to_csv(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a table of gene copy numbers from GFF files.')
    parser.add_argument('output_path', type=str, help='Output path for the gene copy number CSV file.')
    args = parser.parse_args()
    main(args.output_path)




