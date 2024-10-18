from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service
from src.Biocode.utils.utils import tuple_to_sequence_list

@Service
class LinearRepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "linear_repeats_whole_chromosomes"
        self.columns = ["repeats_id", "whole_chromosomes_id", "size", "count"]
        self.pk_column = "id"

    def extract_linear_repeats_by_refseq_accession_number(self, refseq_accession_number: str, k_range: tuple=None):
        return self._extract_repeats_by_method_and_by_refseq_accession_number(refseq_accession_number,
                                                                              method_to_find_it='Linear',
                                                                              size_list=tuple_to_sequence_list(k_range))

    def extract_linear_in_genes_repeats_by_refseq_accession_number(self, refseq_accession_number: str, k_range: tuple=None):
        return self._extract_repeats_by_method_and_by_refseq_accession_number(refseq_accession_number,
                                                                            method_to_find_it='Linear in genes',
                                                                            size_list=tuple_to_sequence_list(k_range))


    def _extract_repeats_by_method_and_by_refseq_accession_number(self, refseq_accession_number: str,
                                                                  method_to_find_it: str, size_list: list):
        query = f"SELECT r.id AS repeats_id, r.name, r.method_to_find_it, lrwc.size, wc.refseq_accession_number, lrwc.count AS count FROM repeats r JOIN " \
                f"linear_repeats_whole_chromosomes lrwc ON r.id = lrwc.repeats_id " \
                f"JOIN whole_chromosomes wc ON lrwc.whole_chromosomes_id = wc.id " \
                f"WHERE r.method_to_find_it = '{method_to_find_it}' AND wc.refseq_accession_number = '{refseq_accession_number}' "

        if size_list:
            query += f"AND lrwc.size IN ({self._list_to_sql_list(size_list)});"
        else:
            query += ";"
        return self.extract_with_custom_query(query)

    def extract_repeats_names_by_size_and_by_refseq_accession_number(self, size: int, refseq_accession_number: str):
        return self._extract_repeats_names_by_size_by_method_and_by_refseq_accession_number(size, 'Linear',
                                                                                            refseq_accession_number)

    def _extract_repeats_names_by_size_by_method_and_by_refseq_accession_number(self, size: int, method_to_find_it: str,
                                                                                refseq_accession_number: str):
        query = f"SELECT r.name, r.id AS repeats_id, lrwc.size, r.method_to_find_it FROM repeats r " \
                f"JOIN linear_repeats_whole_chromosomes lrwc " \
                f"ON r.id = lrwc.repeats_id JOIN whole_chromosomes wc ON lrwc.whole_chromosomes_id = wc.id " \
                f"WHERE lrwc.SIZE = {size} AND r.method_to_find_it = '{method_to_find_it}' " \
                f"AND wc.refseq_accession_number = '{refseq_accession_number}';"
        return self.extract_with_custom_query(query)