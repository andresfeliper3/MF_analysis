import pandas as pd

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
