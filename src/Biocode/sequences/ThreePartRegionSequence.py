from src.Biocode.sequences.RegionSequence import RegionSequence


class ThreePartRegionSequence(RegionSequence):
    def __init__(self, sequence: str = None, sequence_data: dict = None, name: str = None,
                 refseq_accession_number: str = None, organism_name: str = None):
        super().__init__(sequence, sequence_data, 3, name, refseq_accession_number, organism_name)
        super().calculate_regions()

    def set_regions_number(self, regions_number):
        if regions_number != self.regions_number:
            raise Exception("It is not possible to change the regions_1 number of a ThreePartRegionSequence instance")
