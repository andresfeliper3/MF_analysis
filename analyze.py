import os
import sys

project_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..'))  # Get the absolute path of the parent directory
sys.path.append(project_root)  # Add the project root to the Python path

from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from local_scripts.load import data, ORGANISM_NAME, GCF, AMOUNT_CHROMOSOMES, REGIONS_NUMBER
from utils.logger import logger
from utils.timer import timer

from src.Biocode.services.OrganismsService import OrganismsService


@timer
def load_organism(organism_name, gcf, amount_chromosomes):
    print("Loading organism")
    DBConnectionManager.start()
    organism_service = OrganismsService()
    organism_service.insert(record=(organism_name, gcf, amount_chromosomes))
    DBConnectionManager.close()


@timer
def whole_MFA(organism_name, gcf, data):
    DBConnectionManager.start()
    genome_manager = GenomeManager(genome_data=[data[0]], organism_name=organism_name)
    genome_manager.calculate_multifractal_analysis_values()
    genome_manager.save_to_db(GCF=gcf)
    # genome_manager.generate_df_results()

    print(genome_manager.get_mfa_results())
    DBConnectionManager.close()


@timer
def regions_MFA(organism_name, gcf, data, regions_number):
    DBConnectionManager.start()
    region_genome_manager = RegionGenomeManager(genome_data=data, organism_name=organism_name,
                                                regions_number=regions_number)
    region_genome_manager.calculate_multifractal_analysis_values()
    region_genome_manager.save_to_db(GCF=gcf)

    print(region_genome_manager.get_mfa_results())
    DBConnectionManager.close()


load_organism(ORGANISM_NAME, GCF, AMOUNT_CHROMOSOMES)
whole_MFA(ORGANISM_NAME, GCF, data)
