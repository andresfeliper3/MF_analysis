from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class RepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "repeats_whole_chromosomes"
        self.columns = ["repeats_id", "whole_chromosomes_id", "start_position", "end_position", "size"]
        self.pk_column = "id"