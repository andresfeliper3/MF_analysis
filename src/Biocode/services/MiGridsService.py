from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class MiGridsService(AbstractService):
    def __init__(self):
        self.table_name = "mi_grids"
        self.columns = ["mi_grid", "chromosome_id", "epsilon_size"]
        self.pk_column = "id"