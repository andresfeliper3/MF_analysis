from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class KeggCategoriesService(AbstractService):

    def __init__(self):
        self.table_name = "kegg_categories"
        self.columns = ["category"]
        self.pk_column = "id"

    def extract_id_by_category_name(self, category: str) -> int | None:
        result = self.extract_by_field(column='category', value=category)
        if not result.empty:
            return int(result.iloc[0]['id'])
        else:
            return None
