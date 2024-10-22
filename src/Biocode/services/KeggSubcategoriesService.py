from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class KeggSubcategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "kegg_subcategories"
        self.columns = ["subcategory", "category_id"]
        self.pk_column = "id"

    def extract_id_by_subcategory_name(self, subcategory: str) -> int | None:
        result = self.extract_by_field(column='subcategory', value=subcategory)
        if not result.empty:
            return int(result.iloc[0]['id'])
        else:
            return None