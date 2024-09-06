from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.service_decorator import Service

@Service
class WholeMiGridsService(AbstractService):
    def __init__(self):
        self.table_name = "whole_mi_grids"
        self.columns = ["mi_grid", "whole_chromosome_id", "epsilon_size"]
        self.pk_column = "id"

    def extract_mi_grid_by_chromosome_id(self, chromosome_id: int):
        result = self.extract_by_field(column="whole_chromosome_id", value=int(chromosome_id))
        if result is None or result.empty:
            return None
        else:
            return result.loc[0, 'mi_grid']

    def get_chromosome_id_if_mi_grid_exists_by_refseq_accession_number(self, refseq_accession_number):
        query = f"SELECT wc.id AS id FROM whole_mi_grids JOIN whole_chromosomes wc on " \
                f"whole_mi_grids.whole_chromosome_id = wc.id " \
                f"WHERE refseq_accession_number='{refseq_accession_number}';"
        result = self.extract_with_custom_query(query)
        if result is None or result.empty:
            return None
        else:
            return result.loc[0, 'id']
