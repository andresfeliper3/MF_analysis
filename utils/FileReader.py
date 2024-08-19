import pandas as pd
from typing import List

class FileReader:

    @staticmethod
    def read_RM_results_file(path: str):
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
            'refseq_accession_number', 'query_begin', 'query_end', 'query_left', 'strand', 'name',
            'class_family', 'repeat_begin', 'repeat_end', 'repeat_left', 'ID', 'add'
        ])

        df['query_end'] = df['query_end'].astype('int')
        df['query_begin'] = df['query_begin'].astype('int')
        df["repeat_length"] = df["query_end"] - df["query_begin"] + 1
        return df

    @staticmethod
    def divide_genome_df_rows_by_chromosome(df: pd.DataFrame) -> List[pd.DataFrame]:
        df = df[df['sequence'].str.startswith("NC")]
        df_list = []
        start_idx = 0

        for i in range(1, len(df)):
            # If chromosome changes in df
            if df.loc[i, 'sequence'] != df.loc[i - 1, 'sequence']:
                # Append the slice of the DataFrame to the list with index reset
                df_list.append(df.iloc[start_idx:i].reset_index(drop=True))
                start_idx = i

        # Append the last segment with index reset
        df_list.append(df.iloc[start_idx:].reset_index(drop=True))
        return df_list