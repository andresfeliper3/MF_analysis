from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GtfGenesKeggCategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "gtf_genes_kegg_categories"
        self.columns = ["gtf_genes_id", "kegg_categories_id"]
        self.pk_column = "id"

    def extract_count_of_repeats_per_category(self, refseq_accession_number: str):
        query = f"SELECT r.id, r.name, kc.category, gcr.count FROM repeats r " \
                f"JOIN genes_containing_repeats gcr ON r.id = gcr.repeats_id " \
                f"JOIN gtf_genes gg ON gcr.gtf_genes_id=gg.id " \
                f"JOIN gtf_genes_kegg_categories ggkc ON gg.id=ggkc.gtf_genes_id " \
                f"JOIN kegg_categories kc ON kc.id = ggkc.kegg_categories_id " \
                f"JOIN whole_chromosomes wc ON wc.id = gg.whole_chromosomes_id " \
                f"WHERE r.method_to_find_it = 'Linear in genes' " \
                f"AND wc.refseq_accession_number = '{refseq_accession_number}' " \
                f"GROUP BY r.name, kc.category;"
        return self.extract_with_custom_query(query)
