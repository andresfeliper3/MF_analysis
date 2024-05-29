from src.Biocode.services.AbstractService import AbstractService


class RegionChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "region_chromosomes"
        self.columns = ["name", "organism_id", "cover_percentage", "cover", "regions_total", "region_number", "size"]
        self.pk_column = "id"