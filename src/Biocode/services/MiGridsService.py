from src.Biocode.services.AbstractService import AbstractService


class MiGridsService(AbstractService):
    def __init__(self):
        self.table_name = "mi_grids"
        self.columns = ["mi_grid", "chromosome_id", "epsilon_size"]
        self.pk_column = "id"