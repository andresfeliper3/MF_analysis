from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service


@Service
class GenesContainingRepeatsService(AbstractService):

    def __init__(self):
        self.table_name = "genes_containing_repeats"
        self.columns = ["gtf_genes_id", "repeats_id", "count"]
        self.pk_column = "id"