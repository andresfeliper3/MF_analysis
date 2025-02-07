from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service
from src.Biocode.utils.utils import tuple_to_sequence_list


@Service
class LinearRepeatsRegionChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "linear_repeats_region_chromosomes"
        self.columns = ["repeats_id", "region_chromosomes_id", "size", "count"]
        self.pk_column = "id"

    def extract_linear_repeats_by_refseq_accession_number(self, refseq_accession_number: str, k_range: tuple = None):
        return self._extract_repeats_by_method_and_by_refseq_accession_number(refseq_accession_number,
                                                                              method_to_find_it='Linear',
                                                                              size_list=tuple_to_sequence_list(k_range))

    def extract_linear_in_genes_repeats_by_refseq_accession_number(self, refseq_accession_number: str,
                                                                   k_range: tuple = None):
        return self._extract_repeats_by_method_and_by_refseq_accession_number(refseq_accession_number,
                                                                              method_to_find_it='Linear in genes',
                                                                              size_list=tuple_to_sequence_list(k_range))

    def _extract_repeats_by_method_and_by_refseq_accession_number(self, refseq_accession_number: str,
                                                                  method_to_find_it: str, size_list: list):
        query = f"SELECT r.id AS repeats_id, r.name, r.method_to_find_it, lrrc.size, " \
                f"wc.refseq_accession_number, rc.region_number, lrrc.count, rc.size AS window_length " \
                f"FROM repeats r JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id " \
                f"JOIN region_chromosomes rc ON lrrc.region_chromosomes_id = rc.id JOIN whole_chromosomes wc " \
                f"ON rc.whole_chromosome_id = wc.id " \
                f"WHERE r.method_to_find_it = '{method_to_find_it}' AND wc.refseq_accession_number = '{refseq_accession_number}' "
        if size_list:
            query += f"AND lrrc.size IN ({self._list_to_sql_list(size_list)});"
        else:
            query += ";"
        return self.extract_with_custom_query(query)

    def extract_count_of_specific_repeat_by_refseq_accession_number(self, repeat: str, refseq_accession_number: str):
        return self._extract_count_of_specific_repeat_by_method_and_by_refseq_accession_number(repeat, 'Linear',
                                                                                               refseq_accession_number)

    def _extract_count_of_specific_repeat_by_method_and_by_refseq_accession_number(self, repeat: str,
                                                                                   method_to_find_it: str,
                                                                                   refseq_accession_number: str):
        query = f"SELECT r.id AS repeats_id , r.name, r.method_to_find_it, lrrc.size, " \
                f"wc.refseq_accession_number, rc.region_number, lrrc.count, rc.SIZE AS window_length " \
                f"FROM repeats r JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id " \
                f"JOIN region_chromosomes rc ON lrrc.region_chromosomes_id = rc.id " \
                f"JOIN whole_chromosomes wc ON rc.whole_chromosome_id = wc.id " \
                f"WHERE r.method_to_find_it = '{method_to_find_it}' AND wc.refseq_accession_number = '{refseq_accession_number}' " \
                f"AND r.name='{repeat}' " \
                f"ORDER BY rc.region_number;"
        return self.extract_with_custom_query(query)

    def extract_count_of_repeats_by_size_and_chromosome(self, GCF: str, size_list: list):
        query = f"SELECT rc.whole_chromosome_id, SUM(lrrc.count) AS counts_summary FROM repeats r " \
                f"JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id " \
                f"JOIN  region_chromosomes rc ON lrrc.region_chromosomes_id = r.id " \
                f"JOIN organisms o ON o.id = rc.organism_id " \
                f"WHERE r.method_to_find_it = 'Linear' " \
                f"AND o.GCF = '{GCF}' " \
                f"AND lrrc.size IN ({self._list_to_sql_list(size_list)}) " \
                f"GROUP BY rc.whole_chromosome_id;"
        return self.extract_with_custom_query(query)

    def extract_count_of_repeats_by_region(self, refseq_accession_number: str):
        query = f"SELECT rc.id AS region_id , rc.name , SUM(lrrc.count) AS count_sum " \
                f"FROM repeats r JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id " \
                f"JOIN region_chromosomes rc ON lrrc.region_chromosomes_id = rc.id JOIN whole_chromosomes wc ON rc.whole_chromosome_id = wc.id " \
                f"WHERE r.method_to_find_it = 'Linear' AND wc.refseq_accession_number = '{refseq_accession_number}' " \
                f"GROUP BY rc.name ORDER BY rc.id;"
        return self.extract_with_custom_query(query)

    def extract_count_of_single_repeat_by_region(self, refseq_accession_number: str, repeat: str):
        query = f"SELECT rc.id AS region_id , rc.name , SUM(lrrc.count) AS count_sum " \
                f"FROM repeats r JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id " \
                f"JOIN region_chromosomes rc ON lrrc.region_chromosomes_id = rc.id JOIN whole_chromosomes wc ON rc.whole_chromosome_id = wc.id " \
                f"WHERE r.method_to_find_it = 'Linear' AND r.name = '{repeat}' AND wc.refseq_accession_number = '{refseq_accession_number}' " \
                f"GROUP BY rc.name ORDER BY rc.id;"
        return self.extract_with_custom_query(query)

    def extract_most_common_kmer(self, refseq_accession_number_list: list, kmer_size: int) -> str:
        query = (
            f"SELECT r.name, SUM(lrrc.count) AS count "
            f"FROM repeats r "
            f"JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id "
            f"JOIN region_chromosomes rc ON lrrc.region_chromosomes_id = rc.id "
            f"JOIN whole_chromosomes wc ON rc.whole_chromosome_id = wc.id "
            f"WHERE wc.refseq_accession_number IN ({self._list_to_sql_list_with_quotes(refseq_accession_number_list)}) "
            f"AND lrrc.size = {kmer_size} "
            f"GROUP BY r.name "
            f"ORDER BY count DESC "
            f"LIMIT 1;"
        )
        return self.extract_with_custom_query(query).iloc[0]['name']
