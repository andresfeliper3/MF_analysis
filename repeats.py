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
    df = FileReader.read_RM_results_file(path)
    repeats_service = RepeatsService()
    rm_repeats_whole_chromosomes_service= RMRepeatsWholeChromosomesService()
    whole_chromosomes_service = WholeChromosomesService()

    for _, row in df.iterrows():
        repeat_id = repeats_service.insert(record=(row['name'], row['class_family'], 'RM'))
        whole_chromosome_id = whole_chromosomes_service.extract_id_by_refseq_accession_number(row['refseq_accession_number'])
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

@DBConnection
@Timer
def load_RM_repeats_from_folder(path):
    apply_function_to_files_in_folder(path, load_RM_repeats_from_file)