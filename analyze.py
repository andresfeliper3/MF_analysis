from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.OrganismsService import OrganismsService

from load import loader
from utils.timer import timer
from utils.logger import logger


@timer
def load_organism(organism_name, gcf, amount_chromosomes):
    logger.info("Loading organism")
    DBConnectionManager.start()
    organism_service = OrganismsService()
    organism_service.insert(record=(organism_name, gcf, amount_chromosomes))
    DBConnectionManager.close()


@timer
def whole_MFA(organism_name, gcf, data):
    DBConnectionManager.start()
    genome_manager = GenomeManager(genome_data=[data[0]], organism_name=organism_name)
    genome_manager.calculate_multifractal_analysis_values()
    #genome_manager.save_to_db(GCF=gcf)
    # genome_manager.generate_df_results()

    logger.info(genome_manager.get_mfa_results())
    DBConnectionManager.close()


@timer
def regions_MFA(organism_name, gcf, data, regions_number):
    DBConnectionManager.start()
    region_genome_manager = RegionGenomeManager(genome_data=data, organism_name=organism_name,
                                                regions_number=regions_number)
    region_genome_manager.calculate_multifractal_analysis_values()
    region_genome_manager.save_to_db(GCF=gcf)

    logger.info(region_genome_manager.get_mfa_results())
    DBConnectionManager.close()
