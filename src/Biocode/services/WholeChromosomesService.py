from src.Biocode.services.AbstractService import AbstractService


class WholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "whole_chromosomes"
        self.columns = ["name", "organism_id", "cover_percentage", "cover"]
        self.pk_column = "id"