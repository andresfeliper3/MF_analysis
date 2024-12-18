from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GtfGenesKeggSubcategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "gtf_genes_kegg_subcategories"
        self.columns = ["gtf_genes_id", "kegg_subcategories_id"]
        self.pk_column = "id"

    def extract_count_of_repeats_per_subcategory_by_size_and_chromosome(self, refseq_accession_number: str, size: int):
        query = f"SELECT r.id, r.name, ks.subcategory,kc.category, SUM(gcr.count) count FROM repeats r " \
                f"JOIN genes_containing_repeats gcr ON r.id = gcr.repeats_id " \
                f"JOIN gtf_genes gg ON gcr.gtf_genes_id=gg.id " \
                f"JOIN gtf_genes_kegg_subcategories ggks ON gg.id=ggks.gtf_genes_id " \
                f"JOIN kegg_subcategories ks ON ks.id = ggks.kegg_subcategories_id " \
                f"JOIN whole_chromosomes wc ON wc.id = gg.whole_chromosomes_id " \
                f"JOIN kegg_categories kc ON ks.category_id = kc.id " \
                f"WHERE r.method_to_find_it = 'Linear in genes' " \
                f"AND wc.refseq_accession_number = '{refseq_accession_number}' " \
                f"AND LENGTH(r.name) = {size} GROUP BY r.name, ks.subcategory;"
        return self.extract_with_custom_query(query)

    def extract_count_of_repeats_per_subcategory_by_size_and_genome(self, GCF: str, size: int):
        query = f"SELECT r.id, r.name, ks.subcategory, kc.category, SUM(gcr.count) count FROM repeats r " \
                f"JOIN genes_containing_repeats gcr ON r.id = gcr.repeats_id " \
                f"JOIN gtf_genes gg ON gcr.gtf_genes_id=gg.id " \
                f"JOIN gtf_genes_kegg_subcategories ggks ON gg.id=ggks.gtf_genes_id " \
                f"JOIN kegg_subcategories ks ON ks.id = ggks.kegg_subcategories_id " \
                f"JOIN whole_chromosomes wc ON wc.id = gg.whole_chromosomes_id " \
                f"JOIN organisms o ON o.id = wc.organism_id " \
                f"JOIN kegg_categories kc ON ks.category_id = kc.id " \
                f"WHERE r.method_to_find_it = 'Linear in genes' " \
                f"AND o.GCF = '{GCF}' " \
                f"AND LENGTH(r.name) = {size} " \
                f"GROUP BY r.name, ks.subcategory;"
        return self.extract_with_custom_query(query)
