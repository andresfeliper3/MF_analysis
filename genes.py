from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.GtfGenesService import GtfGenesService

import pandas as pd
import sys

from utils.decorators import Timer, DBConnection
from utils.logger import logger
from utils.FileReader import FileReader


@DBConnection
@Timer
def load_genes_from_file(path):
    logger.info(f"Loading genes data from {path}")
    file_df = FileReader.read_gtf_file(path)

    file_df = file_df[file_df['feature'] == 'gene']
    file_df = file_df.reset_index(drop=True)
    
    df_list = FileReader.divide_genome_df_rows_by_chromosome(file_df)
    whole_chromosomes_service = WholeChromosomesService()
    gtf_genes_service = GtfGenesService()

    for df in df_list:
        for _, row in df.iterrows():
            try:
                whole_chromosome_id = whole_chromosomes_service.extract_id_by_refseq_accession_number(
                    row['refseq_accession_number'])
            except Exception as e:
                logger.error(
                    f"Failed to extract whole chromosome ID for refseq_accession_number {row['refseq_accession_number']}: {e}")
                break

            attributes_dict = FileReader.divide_gtf_attributes(row['attributes'])
            record = (
                whole_chromosome_id,
                row['source'],
                row['feature'],
                row['start_position'],
                row['end_position'],
                row['length'],
                row['score'],
                row['strand'],
                row['frame'],
                attributes_dict.get('gene_id', ""),
                attributes_dict.get('gene', ""),
                attributes_dict.get('gene_biotype', "")
            )
            gtf_genes_service.insert(record=record)

