from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GenesContainingRepeatsService(AbstractService):

    def __init__(self):
        self.table_name = "genes_containing_repeats"
        self.columns = ["gtf_genes_id", "repeats_id", "count"]
        self.pk_column = "id"

    def extract_distinct_gene_id_gtf_by_size_and_refseq_accession_number(self, size: int, refseq_accession_number: str):
        query = f"SELECT DISTINCT gg.id, gg.gene_id_gtf FROM genes_containing_repeats gcr JOIN "\
                f"gtf_genes gg ON gcr.gtf_genes_id = gg.id JOIN repeats r ON gcr.repeats_id = r.id "\
                f"JOIN whole_chromosomes wc ON lrwc.whole_chromosomes_id = wc.id "\
                f"JOIN linear_repeats_whole_chromosomes lrwc ON lrwc.repeats_id = r.id "\
                f"WHERE r.method_to_find_it='Linear in genes' AND wc.refseq_accession_number='{refseq_accession_number}' " \
                f"AND lrwc.SIZE = {size} ORDER BY gcr.count DESC"
        return self.extract_with_custom_query(query)