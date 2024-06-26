from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class RegionChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "region_chromosomes"
        self.columns = ["name", "organism_id", "cover_percentage", "cover", "regions_total", "region_number", "size"]
        self.pk_column = "id"