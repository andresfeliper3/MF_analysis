from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager

from load import loader
from utils.timer import timer
from utils.logger import logger


@timer
def load_organism(organism_name, gcf, amount_chromosomes):
    logger.info(f"Loading organism {organism_name} - {gcf} - {amount_chromosomes} chromosomes")
    DBConnectionManager.start()
    organism_service = OrganismsService()
    organism_service.insert(record=(organism_name, gcf, amount_chromosomes))
    DBConnectionManager.close()


@timer
def whole_MFA_genome(organism_name, gcf, data):
    DBConnectionManager.start()
    genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
    genome_manager.calculate_multifractal_analysis_values(GCF=gcf)
    #genome_manager.save_to_db_after_execution(GCF=gcf)
    #genome_manager.graph_linear_fit()
    # genome_manager.generate_df_results()

    DBConnectionManager.close()

@timer
def whole_MFA_sequence(organism_name, sequence_name, gcf, sequence):
    DBConnectionManager.start()
    sequence_manager = SequenceManager(sequence=sequence, organism_name=organism_name,
                                       sequence_name=sequence_name)
    sequence_manager.calculate_multifractal_analysis_values()
    logger.info(sequence_manager.get_mfa_results())
    DBConnectionManager.close()

@timer
def regions_MFA_genome(organism_name, gcf, data, regions_number):
    DBConnectionManager.start()
    region_genome_manager = RegionGenomeManager(genome_data=data, organism_name=organism_name,
                                                regions_number=regions_number)
    region_genome_manager.calculate_multifractal_analysis_values(GCF=gcf)
    #region_genome_manager.save_to_db_after_execution(GCF=gcf)

    DBConnectionManager.close()

@timer
def regions_MFA_sequence(organism_name, sequence_name, gcf, sequence, regions_number):
    DBConnectionManager.start()
    region_sequence_manager = RegionSequenceManager(sequence=sequence, organism_name=organism_name,
                                                regions_number=regions_number, sequence_name=sequence_name)
    region_sequence_manager.calculate_multifractal_analysis_values()
    #region_sequence_manager.save_to_db(GCF=gcf)

    logger.info(region_sequence_manager.get_mfa_results())
    DBConnectionManager.close()
