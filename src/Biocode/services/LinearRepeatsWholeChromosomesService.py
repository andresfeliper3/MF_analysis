from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

@Service
class LinearRepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "linear_repeats_whole_chromosomes"
        self.columns = ["repeats_id", "whole_chromosomes_id", "size"]
        self.pk_column = "id"