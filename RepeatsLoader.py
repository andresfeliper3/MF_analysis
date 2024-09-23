import pandas as pd
from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from utils.FileReader import FileReader
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.folder import apply_function_to_files_in_folder
from utils.logger import logger


@Inject(repeats_service = RepeatsService,
        rm_repeats_service = RMRepeatsWholeChromosomesService,
        whole_chromosomes_service = WholeChromosomesService,
        file_reader=FileReader)
class RepeatsLoader:

    def __init__(self, repeats_service=None, rm_repeats_service=None, whole_chromosomes_service=None, file_reader=None):
        self.repeats_service = repeats_service
        self.rm_repeats_service = rm_repeats_service
        self.whole_chromosomes_service = whole_chromosomes_service
        self.file_reader = file_reader

    @DBConnection
    @TryExcept
    @Timer
    def load_RM_repeats_from_file(self, path):
        logger.info(f"Loading repeats results from RepeatMasker from {path}")
        df = self.file_reader.read_repeats_results_file(path)
        self._load_repeats_file_to_database(df, 'RM')

    @DBConnection
    @TryExcept
    @Timer
    def load_RM_repeats_from_folder(self, path):
        apply_function_to_files_in_folder(path, self.load_RM_repeats_from_file)

    @DBConnection
    @TryExcept
    @Timer
    def load_genome_repeats_file(self, path):
        logger.info(f"Loading repeats results of a genome from a single file in {path}")
        df = self.file_reader.read_repeats_results_file(path)
        df_list = self.file_reader.divide_genome_df_rows_by_chromosome(df)

        for df in df_list:
            self._load_repeats_file_to_database(df, 'Genome in file')

    def _load_repeats_file_to_database(self, df: pd.DataFrame, method_to_find_it: str):
        logger.info(f"*********** Starting loading the chromosome: {df['refseq_accession_number'][0]} ***********")

        for _, row in df.iterrows():
            repeat_id = self.repeats_service.insert(record=(row['repeat'], row['class_family'], method_to_find_it))
            whole_chromosome_id = self._get_whole_chromosome_id(row['refseq_accession_number'])
            if whole_chromosome_id:
                self._insert_rm_repeat(row, repeat_id, whole_chromosome_id)

        logger.info(f"*********** Completed loading the chromosome: {df['refseq_accession_number'][0]} ***********")


    def _get_whole_chromosome_id(self, refseq_accession_number):
        try:
            return self.whole_chromosomes_service.extract_id_by_refseq_accession_number(refseq_accession_number)
        except Exception as e:
            logger.error(
                f"Failed to extract whole chromosome ID for refseq_accession_number {refseq_accession_number}: {e}")
            return None

    def _insert_rm_repeat(self, row, repeat_id, whole_chromosome_id):
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
        self.rm_repeats_service.insert(record=record)
