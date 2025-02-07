import json
import re
from typing import Any

import pandas as pd


def list_to_str(lst):
    return json.dumps(lst)


def str_to_list(string):
    return json.loads(string)


def tuple_to_sequence_list(range_tuple: tuple):
    if range_tuple is None:
        return None
    return list(range(range_tuple[0], range_tuple[1] + 1))


def remove_region_part(input_string: str) -> str:
    match = re.match(r"(.+)_region_\d+$", input_string)
    if match:
        return match.group(1)
    else:
        raise ValueError("Input string does not match the required pattern.")


def adapt_dataframe_to_window_profiles(df: pd.DataFrame) -> list[dict]:
    window_profiles = []
    unique_region_numbers_df_list = df['region_number'].unique()

    for region_number in unique_region_numbers_df_list:
        same_window_df = df[df['region_number'] == region_number]
        window_profile = {}

        for _, row in same_window_df.iterrows():
            kmer = row['name']
            size = row['size']
            count = row['count']

            kmer_key = f"{size}-mers"

            if kmer_key not in window_profile:
                window_profile[kmer_key] = {}

            window_profile[kmer_key][kmer] = count

        window_profiles.append(window_profile)

    return window_profiles


def adapt_dataframe_to_most_frequent_nplets(df: pd.DataFrame) -> dict[str, list[Any]]:
    most_frequent_nplets = {}
    grouped = df.groupby('size')[['name', 'count']].apply(lambda x: list(zip(x['name'], x['count']))).to_dict()

    for size, nplets in grouped.items():
        key = f"{size}-mers"
        sorted_nplets = sorted(nplets, key=lambda x: x[1], reverse=True)
        most_frequent_nplets[key] = sorted_nplets

    return most_frequent_nplets


def get_size_repeats_from_most_frequent_nplets(most_frequent_nplets: dict[str, list], size: int) -> list:
    size_repeats_list = []
    tuples_list = most_frequent_nplets[f"{size}-mers"]
    for kmer_tuple in tuples_list:
        size_repeats_list.append(kmer_tuple[0])
    return size_repeats_list

