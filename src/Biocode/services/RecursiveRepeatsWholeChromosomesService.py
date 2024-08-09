from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

from typing import List

@Service
class RecursiveRepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "recursive_repeats_whole_chromosomes"
        self.columns = ["repeats_id", "whole_chromosomes_id", "size", "largest_value", "coordinates"]
        self.pk_column = "id"

    def extract_info_by_chromosome(self, refseq_accession_number: str):
        query = f"SELECT repeats.name AS name, rwc.size AS repeat_length, largest_value, coordinates FROM repeats " \
                f"JOIN recursive_repeats_whole_chromosomes rwc on repeats.id = rwc.repeats_id " \
                f"LEFT JOIN whole_chromosomes wc on rwc.whole_chromosomes_id = wc.id " \
                f"WHERE refseq_accession_number='{refseq_accession_number}';"
        return DBConnectionManager.extract_with_custom_query(query)


