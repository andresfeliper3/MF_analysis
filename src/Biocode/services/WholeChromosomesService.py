from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class WholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "whole_chromosomes"
        self.columns = ["name", "organism_id", "cover_percentage", "cover", "size"]
        self.pk_column = "id"