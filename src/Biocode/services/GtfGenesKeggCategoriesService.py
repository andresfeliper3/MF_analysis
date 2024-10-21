from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GtfGenesKeggCategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "gtf_genes_kegg_categories"
        self.columns = ["gtf_genes_id", "kegg_categories_id"]
        self.pk_column = "id"
