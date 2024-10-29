import pandas as pd
from typing import List, Dict


class FileReader:

    @staticmethod
    def read_repeats_results_file(path: str):
        data = []

        with open(path, 'r') as file:
            next(file)  # Skip 3 lines
            next(file)
            next(file)
            for line in file:
                # Split the line by whitespace(s)
                columns = line.strip().split()
                # Append the columns to the data list
                data.append(columns)

        df = pd.DataFrame(data, columns=[
            'sw_score', 'percentage_divergence', 'percentage_deletions', 'percentage_insertions',
            'refseq_accession_number', 'query_begin', 'query_end', 'query_left', 'strand', 'repeat',
            'class_family', 'repeat_begin', 'repeat_end', 'repeat_left', 'ID', 'add'
        ])

        df['query_end'] = df['query_end'].astype('int')
        df['query_begin'] = df['query_begin'].astype('int')
        df["repeat_length"] = df["query_end"] - df["query_begin"] + 1
        return df

    @staticmethod
    def divide_genome_df_rows_by_chromosome(df: pd.DataFrame) -> List[pd.DataFrame]:
        df = df[df['refseq_accession_number'].str.startswith("NC")]
        df_list = []
        start_idx = 0
        for i in range(1, len(df)):
            # If chromosome changes in df
            if df.iloc[i]['refseq_accession_number'] != df.iloc[i - 1]['refseq_accession_number']:
                # Append the slice of the DataFrame to the list with index reset
                df_list.append(df.iloc[start_idx:i].reset_index(drop=True))
                start_idx = i

        # Append the last segment with index reset
        df_list.append(df.iloc[start_idx:].reset_index(drop=True))
        return df_list

    @staticmethod
    def read_gtf_file(path):
        # Define column names for GFF file
        column_names = ["refseq_accession_number", "source", "feature", "start_position", "end_position", "score", "strand", "frame", "attributes"]

        df = pd.read_csv(path, sep='\t', comment='#', header=None, names=column_names)

        """
        # Add a new column "Chromosome" initialized with 1
        df.insert(0, "Chromosome", 1)

        # Iterate over the DataFrame to update "Chromosome" values
        current_chromosome = 1
        for index in range(1, len(df)):
            if df.iloc[index]['refseq_accession_number'] != df.iloc[index - 1]['refseq_accession_number']:
                current_chromosome += 1
            df.at[df.index[index], 'Chromosome'] = current_chromosome
        """
        # Convert 'start' and 'end' columns to integers
        df['start_position'] = df['start_position'].astype(int)
        df['end_position'] = df['end_position'].astype(int)

        # Calculate 'length' column
        df["length"] = df["end_position"] - df["start_position"] + 1

        return df

    @staticmethod
    def divide_gtf_attributes(attributes: str) -> Dict[str, str]:
        keys_of_interest = {"gene_id", "gene", "gene_biotype"}
        split_attributes = attributes.split(';')
        attribute_dict = {}
        # Iterate through each key-value pair
        for attribute in split_attributes:
            if attribute.strip():  # Make sure it's not an empty string
                try:
                    key, value = attribute.strip().split(' ', 1)
                    key = key.strip()
                    value = value.strip('"')
                except ValueError:
                    # Handle the case where split does not return two values
                    key = attribute.strip()
                    value = "true"

                if key in keys_of_interest:
                    attribute_dict[key] = value

        return attribute_dict
