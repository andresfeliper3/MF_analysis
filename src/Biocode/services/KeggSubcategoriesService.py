from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class KeggSubcategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "kegg_subcategories"
        self.columns = ["subcategory", "category_id"]
        self.pk_column = "id"
