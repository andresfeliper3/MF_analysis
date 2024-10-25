from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GtfGenesKeggSubcategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "gtf_genes_kegg_subcategories"
        self.columns = ["gtf_genes_id", "kegg_subcategories_id"]
        self.pk_column = "id"

    def extract_count_of_repeats_per_subcategory(self, refseq_accession_number: str):
        query = f"SELECT r.id, r.name, ks.subcategory, gcr.count FROM repeats r " \
                f"JOIN genes_containing_repeats gcr ON r.id = gcr.repeats_id " \
                f"JOIN gtf_genes gg ON gcr.gtf_genes_id=gg.id " \
                f"JOIN gtf_genes_kegg_subcategories ggks ON gg.id=ggks.gtf_genes_id " \
                f"JOIN kegg_subcategories ks ON ks.id = ggks.kegg_subcategories_id " \
                f"JOIN whole_chromosomes wc ON wc.id = gg.whole_chromosomes_id " \
                f"WHERE r.method_to_find_it = 'Linear in genes' " \
                f"AND wc.refseq_accession_number = '{refseq_accession_number}' " \
                f"GROUP BY r.name, ks.subcategory;"
