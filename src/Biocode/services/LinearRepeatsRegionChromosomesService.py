from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class LinearRepeatsRegionChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "linear_repeats_region_chromosomes"
        self.columns = ["repeats_id", "region_chromosomes_id", "size", "count"]
        self.pk_column = "id"

    def extract_linear_repeats_by_refseq_accession_number(self, refseq_accession_number: str):
        return self._extract_repeats_by_method_and_by_refseq_accession_number(refseq_accession_number,
                                                                              method_to_find_it='Linear')

    def _extract_repeats_by_method_and_by_refseq_accession_number(self, refseq_accession_number: str,
                                                                  method_to_find_it: str):
        query = f"SELECT r.id AS repeats_id, r.name, r.method_to_find_it, lrrc.size, " \
                 f"wc.refseq_accession_number, rc.region_number, lrrc.count " \
                 f"FROM repeats r JOIN linear_repeats_region_chromosomes lrrc ON r.id = lrrc.repeats_id " \
                 f"JOIN region_chromosomes rc ON lrrc.region_chromosomes_id = rc.id JOIN whole_chromosomes wc " \
                 f"ON rc.whole_chromosome_id = wc.id " \
                 f"WHERE r.method_to_find_it = '{method_to_find_it}' AND wc.refseq_accession_number = '{refseq_accession_number}';"
        return self.extract_with_custom_query(query)
