from src.Biocode.services.WholeResultsService import WholeResultsService
from src.Biocode.services.RegionResultsService import RegionResultsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RecursiveRepeatsWholeChromosomesService import RecursiveRepeatsWholeChromosomesService
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.graphs.Graphs import Graphs
from load import loader
from src.Biocode.utils.utils import str_to_list
from utils.decorators import Timer, DBConnection
from utils.logger import logger


@Timer
def load_data_whole(gcf) -> dict:
    DBConnectionManager.start()
    whole_results_service = WholeResultsService()
    df = whole_results_service.extract_results(GCF=gcf)
    DBConnectionManager.close()
    return df.to_dict(orient='records')


@Timer
def graph_whole(dataframe, organism_name, data):
    genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
    desired_keys_ddq = ['DDq', 'sequence_name']
    desired_keys_dq_tauq = ['Dq_values', 'tau_q_values']

    mfa_results = []

    for item in dataframe:
        result_entry = {
            'q_values': list(range(-20, 21)),
            **{key: item[key] for key in desired_keys_ddq}
        }

        result_entry.update({
            key: str_to_list(item[key]) for key in desired_keys_dq_tauq
        })

        mfa_results.append(result_entry)

    cover = [str_to_list(item['cover']) for item in dataframe]
    cover_percentage = [item['cover_percentage'] for item in dataframe]
    degrees_of_multifractality = [item['DDq'] for item in dataframe]

    genome_manager.set_mfa_results(mfa_results)
    genome_manager.set_cover(cover)
    genome_manager.set_cover_percentage(cover_percentage)
    genome_manager.set_degrees_of_multifractality(degrees_of_multifractality)

    genome_manager.generate_df_results()

    genome_manager.graph_degrees_of_multifractality()
    genome_manager.graph_multifractal_analysis_merged()

    genome_manager.graph_coverage()


@Timer
def load_data_regions(gcf) -> dict:
    DBConnectionManager.start()
    region_results_service = RegionResultsService()
    df = region_results_service.extract_results(GCF=gcf)
    DBConnectionManager.close()
    return df.to_dict(orient='records')


@Timer
def graph_regions(dataframe, organism_name, data, regions_number):
    region_genome_manager = RegionGenomeManager(genome_data=data, organism_name=organism_name,
                                                regions_number=regions_number)

    desired_keys_ddq = ['DDq', 'sequence_name']
    desired_keys_dq_tauq = ['Dq_values', 'tau_q_values']

    mfa_results = []

    for item in dataframe:
        result_entry = {
            'q_values': list(range(-20, 21)),
            **{key: item[key] for key in desired_keys_ddq}
        }

        result_entry.update({
            key: str_to_list(item[key]) for key in desired_keys_dq_tauq
        })

        mfa_results.append(result_entry)

    cover = [str_to_list(item['cover']) for item in dataframe]
    cover_percentage = [item['cover_percentage'] for item in dataframe]
    degrees_of_multifractality = [item['DDq'] for item in dataframe]

    region_genome_manager.set_mfa_results(mfa_results)
    region_genome_manager.set_flattened_mfa_results(mfa_results)
    region_genome_manager.set_cover(cover)
    region_genome_manager.set_cover_percentage(cover_percentage)
    region_genome_manager.set_degrees_of_multifractality(degrees_of_multifractality)

    region_genome_manager.generate_df_results()

    region_genome_manager.graph_degrees_of_multifractality()
    region_genome_manager.graph_multifractal_analysis_merged()

    region_genome_manager.graph_coverage()

@DBConnection
@Timer
def graph_rm_results_from_file(path: str, refseq_accession_number:str, partitions:int, regions: int,
                               plot_type:str = "line", save: bool = True, name: str = None):
    DEFAULT_REGIONS = 3
    DEFAULT_PARTITIONS = 300
    DEFAULT_REPEATS_LIMIT = 20
    whole_chromosomes_service = WholeChromosomesService()
    size = whole_chromosomes_service.extract_size_by_refseq_accession_number(refseq_accession_number)
    partitions = int(partitions) if isinstance(partitions, str) else DEFAULT_PARTITIONS
    regions = int(regions) if isinstance(regions, str) else DEFAULT_REGIONS
    Graphs.graph_distribution_of_repeats_merged_from_file(path=path, size=size, partitions=partitions,
                                             legend=True, regions=regions, plot_type=plot_type, save=save,
                                             name=name, refseq_accession_number=refseq_accession_number)


    Graphs.graph_frequency_of_repeats_grouped_from_file(path, col="class_family", filtering=False, n_max=10, save=save,
                                                        name=name, refseq_accession_number=refseq_accession_number)
    Graphs.graph_frequency_of_repeats_grouped_from_file(path, col="name", filtering=False, n_max=10, save=save, name=name,
                                               refseq_accession_number=refseq_accession_number)

    Graphs.graph_distribution_of_repeats_from_file(path, col="class_family", legend=True, plot_type=plot_type,
                                                   limit=20, regions=regions, save=save, name=name,
                                                   refseq_accession_number=refseq_accession_number)
    Graphs.graph_distribution_of_repeats_from_file(path, col="name", legend=True, plot_type=plot_type,
                                                   limit=20, regions=regions, save=save, name=name,
                                                   refseq_accession_number=refseq_accession_number)

    Graphs.graph_distribution_of_repeats_subplots_from_file(path, col="class_family", legend=True,
                                                            limit=DEFAULT_REPEATS_LIMIT, regions=regions, save=save,
                                                            name=name, refseq_accession_number=refseq_accession_number)
    Graphs.graph_distribution_of_repeats_subplots_from_file(path, col="name", legend=True,
                                                            limit=DEFAULT_REPEATS_LIMIT, regions=regions, save=save,
                                                            name=name, refseq_accession_number=refseq_accession_number)


@DBConnection
@Timer
def graph_rm_results_from_database(refseq_accession_number:str, partitions:int, regions: int,
                               plot_type:str = "line", save: bool = True, name: str = None):
    DEFAULT_REGIONS = 3
    DEFAULT_PARTITIONS = 300
    DEFAULT_REPEATS_LIMIT = 20
    whole_chromosomes_service = WholeChromosomesService()
    size = whole_chromosomes_service.extract_size_by_refseq_accession_number(refseq_accession_number)
    partitions = int(partitions) if isinstance(partitions, str) else DEFAULT_PARTITIONS
    regions = int(regions) if isinstance(regions, str) else DEFAULT_REGIONS

    rm_repeats_whole_chromosomes_service = RMRepeatsWholeChromosomesService()

    data = rm_repeats_whole_chromosomes_service.extract_info_by_chromosome(refseq_accession_number)

    Graphs.graph_distribution_of_repeats_merged_from_database(data=data, size=size, partitions=partitions,
                                             legend=True, regions=regions, plot_type=plot_type, save=save,
                                             name=name, refseq_accession_number=refseq_accession_number)

    Graphs.graph_frequency_of_repeats_grouped_from_database(data, col="class_family", filtering=False, n_max=10, save=save,
                                                        name=name, refseq_accession_number=refseq_accession_number)
    Graphs.graph_frequency_of_repeats_grouped_from_database(data, col="name", filtering=False, n_max=10, save=save, name=name,
                                               refseq_accession_number=refseq_accession_number)

    Graphs.graph_distribution_of_repeats_from_database(data, col="class_family", legend=True, plot_type=plot_type,
                                                   limit=20, regions=regions, save=save, name=name,
                                                   refseq_accession_number=refseq_accession_number)
    Graphs.graph_distribution_of_repeats_from_database(data, col="name", legend=True, plot_type=plot_type,
                                                       limit=20, regions=regions, save=save, name=name,
                                                       refseq_accession_number=refseq_accession_number)

    Graphs.graph_distribution_of_repeats_subplots_from_database(data, col="class_family", legend=True,
                                                                limit=DEFAULT_REPEATS_LIMIT, regions=regions, save=save,
                                                                name=name,
                                                                refseq_accession_number=refseq_accession_number)
    Graphs.graph_distribution_of_repeats_subplots_from_database(data, col="name", legend=True,
                                                           limit=DEFAULT_REPEATS_LIMIT, regions=regions, save=save,
                                                           name=name, refseq_accession_number=refseq_accession_number)


@DBConnection
@Timer
def graph_recursive_from_database(refseq_accession_number: str, save: bool, name: str, n_max: int):
    n_max = n_max and int(n_max)
    recursive_repeats_whole_chromosomes_service = RecursiveRepeatsWholeChromosomesService()

    data = recursive_repeats_whole_chromosomes_service.extract_info_by_chromosome(refseq_accession_number)
    Graphs.graph_recursive_repeats_largest_values_from_database(data, col="name", n_max=n_max, save=save, name=name,
                                                                refseq_accession_number=refseq_accession_number)
    Graphs.graph_grouped_by_recursive_repeat_length(data, col="name", save=save, name=name,
                                                    refseq_accession_number=refseq_accession_number)
    Graphs.graph_individual_plots_by_recursive_repeat_length(data, col="name", save=save, name=name,
                                                    refseq_accession_number=refseq_accession_number)