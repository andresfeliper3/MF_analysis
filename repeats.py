from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService

from load import loader
from utils.decorators import Timer, DBConnection
from utils.logger import logger
from utils.FileReader import FileReader
from utils.folder import apply_function_to_files_in_folder

import pandas as pd


@DBConnection
@Timer
def load_RM_repeats_from_file(path):
    logger.info(f"Loading repeats results from RepeatMasker from {path}")
    df = FileReader.read_repeats_results_file(path)
    repeats_service = RepeatsService()
    rm_repeats_whole_chromosomes_service= RMRepeatsWholeChromosomesService()
    whole_chromosomes_service = WholeChromosomesService()
    _load_repeats_file_to_database(df, 'RM', repeats_service, rm_repeats_whole_chromosomes_service,
                                   whole_chromosomes_service)


@DBConnection
@Timer
def load_RM_repeats_from_folder(path):
    apply_function_to_files_in_folder(path, load_RM_repeats_from_file)


@DBConnection
@Timer
def load_genome_repeats_file(path):
    logger.info(f"Loading repeats results of a genome from a single file in {path}")
    df = FileReader.read_repeats_results_file(path)
    df_list = FileReader.divide_genome_df_rows_by_chromosome(df)
    repeats_service = RepeatsService()
    rm_repeats_whole_chromosomes_service = RMRepeatsWholeChromosomesService()
    whole_chromosomes_service = WholeChromosomesService()

    for df in df_list:
        _load_repeats_file_to_database(df, 'Genome in file', repeats_service, rm_repeats_whole_chromosomes_service,
                                       whole_chromosomes_service)


def _load_repeats_file_to_database(df: pd.DataFrame,  method_to_find_it: str, repeats_service,
                                   rm_repeats_whole_chromosomes_service, whole_chromosomes_service):
    logger.info(f"*********** Starting loading the chromosome: {df['refseq_accession_number'][0]} ***********")
    for _, row in df.iterrows():
        repeat_id = repeats_service.insert(record=(row['repeat'], row['class_family'], method_to_find_it))
        whole_chromosome_id = whole_chromosomes_service.extract_id_by_refseq_accession_number(
            row['refseq_accession_number'])
        record = (
            repeat_id,
            whole_chromosome_id,
            row['sw_score'],
            row['percentage_divergence'],
            row['percentage_deletions'],
            row['percentage_insertions'],
            row['query_begin'],
            row['query_end'],
            row['repeat_length'],
            row['query_left'],
            row['strand'],
            row['repeat_begin'],
            row['repeat_end'],
            row['repeat_left']
        )
        rm_repeats_whole_chromosomes_service.insert(record=record)
    logger.info(f"*********** Completed loading the chromosome: {df['refseq_accession_number'][0]} ***********")