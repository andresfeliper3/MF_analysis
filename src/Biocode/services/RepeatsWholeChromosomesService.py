from src.Biocode.services.AbstractService import AbstractService

class RepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "repeats_whole_chromosomes"
        self.columns = ["repeats_id", "whole_chromosomes_id", "start_position", "end_position", "size"]
        self.pk_column = "id"