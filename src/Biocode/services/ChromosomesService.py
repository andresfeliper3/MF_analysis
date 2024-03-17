from src.Biocode.services.AbstractService import AbstractService


class ChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "chromosomes"
        self.columns = ["name", "organism_id", "cover_percentage", "cover"]
        self.pk_column = "id"