from src.Biocode.services.AbstractService import AbstractService
from src.Biocode.services.services_context.all_services import Service

@Service
class RMRepeatsWholeChromosomesService(AbstractService):
    def __init__(self):
        self.table_name = "RM_repeats_whole_chromosomes"
        self.columns = [
            "repeats_id", "whole_chromosomes_id", "sw_score", "percentage_divergence", "percentage_deletions",
            "percentage_insertions", "query_begin", "query_end", "repeat_length", "query_left", "strand",
            "repeat_begin", "repeat_end", "repeat_left"]

        self.pk_column = "id"