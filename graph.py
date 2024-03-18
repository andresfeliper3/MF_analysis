from src.Biocode.services.WholeResultsService import WholeResultsService
from src.Biocode.services.RegionResultsService import RegionResultsService
from src.Biocode.managers.DBConnectionManager import DBConnectionManager
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from load import loader
from src.Biocode.utils.utils import str_to_list
from utils.timer import timer
from utils.logger import logger


@timer
def load_data_whole(gcf) -> dict:
    DBConnectionManager.start()
    whole_results_service = WholeResultsService()
    df = whole_results_service.extract_results(GCF=gcf)
    DBConnectionManager.close()
    return df.to_dict(orient='records')


@timer
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


@timer
def load_data_regions(gcf) -> dict:
    DBConnectionManager.start()
    region_results_service = RegionResultsService()
    df = region_results_service.extract_results(GCF=gcf)
    DBConnectionManager.close()
    return df.to_dict(orient='records')


@timer
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