from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

@Service
class LinearRepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "linear_repeats_whole_chromosomes"
        self.columns = ["repeats_id", "whole_chromosomes_id", "size"]
        self.pk_column = "id"

    def extract_linear_repeats_by_refseq_accession_number(self, refseq_accession_number: str):
       return self._extract_repeats_by_method_and_by_refseq_accession_number(refseq_accession_number,
                                                                             method_to_find_it='Linear')

    def _extract_repeats_by_method_and_by_refseq_accession_number(self, refseq_accession_number: str, method_to_find_it: str):
        query = f"SELECT r.id AS repeats_id, r.name, r.class_family, r.method_to_find_it, lrwc.size, wc.refseq_accession_number  FROM repeats r JOIN " \
                f"linear_repeats_whole_chromosomes lrwc ON r.id = lrwc.repeats_id " \
                f"JOIN whole_chromosomes wc ON lrwc.whole_chromosomes_id = wc.id " \
                f"WHERE r.method_to_find_it = '{method_to_find_it}' AND refseq_accession_number = '{refseq_accession_number}';"
        return self.extract_with_custom_query(query)

