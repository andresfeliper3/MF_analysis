from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager

from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.RecursiveRepeatsWholeChromosomesService import RecursiveRepeatsWholeChromosomesService

from load import loader
from utils.decorators import Timer, DBConnection
from utils.logger import logger


@DBConnection
@Timer
def load_organism(organism_name, gcf, amount_chromosomes):
    logger.info(f"Loading organism {organism_name} - {gcf} - {amount_chromosomes} chromosomes - to database")
    organism_service = OrganismsService()
    organism_service.insert(record=(organism_name, gcf, amount_chromosomes))


@DBConnection
@Timer
def whole_MFA_genome(organism_name, gcf, data, save_to_db):
    genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
    genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db)
    #genome_manager.save_to_db_after_execution(GCF=gcf)
    #genome_manager.graph_linear_fit()
    # genome_manager.generate_df_results()

@DBConnection
@Timer
def find_kmers_recursively_in_genome(organism_name, gcf, data, save_to_db):
    genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
    genome_manager.find_only_kmers_recursively_and_calculate_multifractal_analysis_values(
        GCF=gcf, save_to_db=save_to_db, method_to_find_it="Recursively")


@DBConnection
@Timer
def whole_MFA_sequence(gcf, sequence, save_to_db):
    sequence_manager = SequenceManager(sequence=sequence)
    sequence_manager.calculate_multifractal_analysis_values()

    if save_to_db:
        sequence_manager.save_to_db_during_execution(GCF=gcf)

    logger.info(sequence_manager.get_mfa_results())

@DBConnection
@Timer
def find_kmers_recursively_in_sequence(gcf, sequence, save_to_db):
    sequence_manager = SequenceManager(sequence=sequence)
    sequence_manager.calculate_multifractal_analysis_values()
    kmers_list = sequence_manager.find_only_kmers_recursively()
    if save_to_db:
        sequence_manager.save_repeats_found_recursively_to_db(
            kmers_list=kmers_list, GCF=gcf, method_to_find_it="Recursively")


@DBConnection
@Timer
def regions_MFA_genome(organism_name, gcf, data, regions_number, save_to_db):
    region_genome_manager = RegionGenomeManager(genome_data=data, regions_number=regions_number,
                                                organism_name=organism_name)
    region_genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db)
    #region_genome_manager.save_to_db_after_execution(GCF=gcf)


@DBConnection
@Timer
def regions_MFA_sequence(gcf, sequence, regions_number, save_to_db):

    region_sequence_manager = RegionSequenceManager(sequence=sequence, regions_number=regions_number)
    region_sequence_manager.calculate_multifractal_analysis_values()

    if save_to_db:
        region_sequence_manager.save_to_db_during_execution(GCF=gcf)

    logger.info(region_sequence_manager.get_mfa_results())
