from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class RegionChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "region_chromosomes"
        self.columns = ["name", "refseq_accession_number", "organism_id", "cover_percentage", "cover", "regions_total", "region_number",
                        "size", "whole_chromosome_id"]
        self.pk_column = "id"