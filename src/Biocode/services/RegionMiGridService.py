from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

@Service
class RegionMiGridsService(AbstractService):
    def __init__(self):
        self.table_name = "region_mi_grids"
        self.columns = ["mi_grid", "region_chromosome_id", "epsilon_size"]
        self.pk_column = "id"

    def extract_mi_grid_by_chromosome_id(self, id: int):
        return self.extract_by_field(column="region_chromosome_id", value=id).loc[0, 'mi_grid']
