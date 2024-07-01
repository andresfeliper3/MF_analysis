from src.Biocode.services.AbstractService import AbstractService

class RepeatsService(AbstractService):

    def __init__(self):
        self.table_name = "repeats"
        self.columns = ["name", "class_family", "method_to_find_it"]
        self.pk_column = "id"