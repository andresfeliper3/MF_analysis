from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class KeggCategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "kegg_categories"
        self.columns = ["category"]
        self.pk_column = "id"
