import yaml

from load import loader
from src.Biocode.graphs.Graphs import Graphs
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RecursiveRepeatsWholeChromosomesService import RecursiveRepeatsWholeChromosomesService
from src.Biocode.services.RegionResultsService import RegionResultsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.WholeResultsService import WholeResultsService
from src.Biocode.utils.utils import str_to_list
from utils.FileReader import FileReader
from utils.decorators import Timer, DBConnection, TryExcept
from utils.folder import apply_function_to_files_in_folder
from utils.logger import logger


@Inject(region_genome_manager = RegionGenomeManager,
        whole_results_service = WholeResultsService,
        region_results_service = RegionResultsService,
        whole_chromosomes_service = WholeChromosomesService,
        recursive_repeats_service = RecursiveRepeatsWholeChromosomesService,
        organisms_service = OrganismsService,
        gtf_genes_service = GtfGenesService)
class Grapher:

    DEFAULT_REGIONS = 3
    DEFAULT_PARTITIONS = 300
    DEFAULT_REPEATS_LIMIT = 20

    def __init__(self,
                 whole_results_service: WholeResultsService, region_results_service: RegionResultsService,
                 whole_chromosomes_service: WholeChromosomesService, recursive_repeats_service: RecursiveRepeatsWholeChromosomesService,
                 organisms_service: OrganismsService, gtf_genes_service: GtfGenesService):
        self.whole_results_service = whole_results_service
        self.region_results_service = region_results_service
        self.whole_chromosomes_service = whole_chromosomes_service
        self.recursive_repeats_service = recursive_repeats_service
        self.organisms_service = organisms_service
        self.gtf_genes_service = gtf_genes_service
        self.graphs_config_path = "config/graphs_config.yaml"
        self.organism = ""

    def load_config(self):
        with open(self.graphs_config_path, 'r') as file:
            return yaml.safe_load(file)

    @DBConnection
    @TryExcept
    @Timer
    def graph_command(self, args):
        if args.name:
            self.organism = args.name
            loader.set_organism(self.organism)
            self._validate_mode_graphing(args)
        else:
            raise Exception("Please provide either -id or -name.")

    def _validate_mode_graphing(self, args):
        if args.mode:
            if args.mode == 'whole':
                dic = self._load_data_whole(gcf=loader.get_gcf())
                self._graph_whole(dataframe=dic, organism_name=loader.get_organism_name(), data=loader.get_data())
            elif args.mode == 'regions':
                dic_list = self._load_data_regions(gcf=loader.get_gcf())
                self._graph_regions(dataframe=dic_list, organism_name=loader.get_organism_name(), data=loader.get_data(),
                                    regions_number=loader.get_regions_number())
            else:
                raise Exception("Enter a valid mode (whole or regions)")
        else:
            raise Exception("Enter a valid mode (whole or regions)")


    def _load_data_whole(self, gcf) -> dict:
        df = self.whole_results_service.extract_results(GCF=gcf)
        return df.to_dict(orient='records')

    def _graph_whole(self, dataframe, organism_name, data):
        config = self.load_config()
        graphs_config = config.get('MFA', {})

        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        mfa_results = self._prepare_mfa_results(dataframe)

        cover, cover_percentage, degrees_of_multifractality = self._extract_cover_info(dataframe)

        genome_manager.set_mfa_results(mfa_results)
        genome_manager.set_cover(cover)
        genome_manager.set_cover_percentage(cover_percentage)
        genome_manager.set_degrees_of_multifractality(degrees_of_multifractality)

        genome_manager.generate_df_results()

        self._graph_MFA_options(graphs_config, genome_manager)

    def _load_data_regions(self, gcf) -> list[dict]:
        region_results_service = RegionResultsService()
        df = region_results_service.extract_results(GCF=gcf)
        return df.to_dict(orient='records')

    def _graph_regions(self, dataframe, organism_name, data, regions_number):
        config = self.load_config()
        graphs_config = config.get('MFA', {})

        region_genome_manager = RegionGenomeManager(
            genome_data=data, organism_name=organism_name, regions_number=regions_number
        )
        mfa_results = self._prepare_mfa_results(dataframe)

        cover, cover_percentage, degrees_of_multifractality = self._extract_cover_info(dataframe)

        region_genome_manager.set_mfa_results(mfa_results)
        region_genome_manager.set_flattened_mfa_results(mfa_results)
        region_genome_manager.set_cover(cover)
        region_genome_manager.set_cover_percentage(cover_percentage)
        region_genome_manager.set_degrees_of_multifractality(degrees_of_multifractality)

        region_genome_manager.generate_df_results()

        self._graph_MFA_options(graphs_config, region_genome_manager)

    def _prepare_mfa_results(self, dataframe):
        """Extract MFA results and return a list of dictionaries"""
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
        return mfa_results

    def _extract_cover_info(self, dataframe):
        cover = [str_to_list(item['cover']) for item in dataframe]
        cover_percentage = [item['cover_percentage'] for item in dataframe]
        degrees_of_multifractality = [item['DDq'] for item in dataframe]
        return cover, cover_percentage, degrees_of_multifractality

    def _graph_MFA_options(self, graphs_config, manager):
        if graphs_config.get('degrees_of_multifractality', False):
            manager.graph_degrees_of_multifractality()

        if graphs_config.get('multifractal_analysis_merged', False):
            manager.graph_multifractal_analysis_merged()

        if graphs_config.get('coverage', False):
            manager.graph_coverage()

    @DBConnection
    @TryExcept
    @Timer
    def graph_rm_results_from_file(self, path: str, partitions: int, regions: int, plot_type: str = "line",
                                   save: bool = True, dir: str = None):
        config = self.load_config()
        graphs_config = config.get('repeats', {})

        df = FileReader.read_repeats_results_file(path)
        refseq_accession_number = df['refseq_accession_number'][0]
        filename, size = self._get_filename_and_size(refseq_accession_number)

        partitions, regions = self._set_default_values(partitions, regions)

        self._graph_rm_options(graphs_config, df, size, partitions, regions, plot_type, save, dir, filename)

    def _get_filename_and_size(self, refseq_accession_number):
        return self.whole_chromosomes_service.extract_filename_and_size_by_refseq_accession_number(refseq_accession_number)

    def _set_default_values(self, partitions, regions):
        partitions = int(partitions) if isinstance(partitions, str) else self.DEFAULT_PARTITIONS
        regions = int(regions) if isinstance(regions, str) else self.DEFAULT_REGIONS
        return partitions, regions

    def _graph_rm_options(self, graphs_config, df, size, partitions, regions, plot_type, save, dir, filename):
        if graphs_config.get('distribution_of_repeats_merged', False):
            Graphs.graph_distribution_of_repeats_merged(df, size, partitions, regions, plot_type, save, dir, filename)

        if graphs_config.get('frequency_of_repeats_class_family', False):
            Graphs.graph_frequency_of_repeats_grouped(df, col="class_family", n_max=10, save=save, name=dir,
                                                      filename=filename)

        if graphs_config.get('frequency_of_repeats_repeat', False):
            Graphs.graph_frequency_of_repeats_grouped(
                data=df, col="repeat", filtering=False, n_max=10,
                save=save, name=dir, filename=filename)

        if graphs_config.get('distribution_of_repeats_class_family', False):
            Graphs.graph_distribution_of_repeats(
                df=df, col="class_family", legend=True, plot_type=plot_type,
                limit=DEFAULT_REPEATS_LIMIT, regions=regions, save=save,
                name=dir, filename=filename)

        if graphs_config.get('distribution_of_repeats_repeat', False):
            Graphs.graph_distribution_of_repeats(
                df=df, col="repeat", legend=True, plot_type=plot_type,
                limit=DEFAULT_REPEATS_LIMIT, regions=regions, save=save,
                name=dir, filename=filename)

        if graphs_config.get('distribution_of_repeats_subplots_class_family', False):
            Graphs.graph_distribution_of_repeats_subplots(
                df=df, col="class_family", legend=True, limit=DEFAULT_REPEATS_LIMIT,
                regions=regions, save=save, name=dir, filename=filename)

        if graphs_config.get('distribution_of_repeats_subplots_repeat', False):
            Graphs.graph_distribution_of_repeats_subplots(
                df=df, col="repeat", legend=True, limit=DEFAULT_REPEATS_LIMIT,
                regions=regions, save=save, name=dir, filename=filename)

    @DBConnection
    @TryExcept
    @Timer
    def graph_rm_results_from_database(self, refseq_accession_number: str, partitions: int, regions: int,
                                       plot_type: str = "line", save: bool = True, dir: str = None):
        config = self.load_config()
        graphs_config = config.get('repeats', {})

        filename, size = self._get_filename_and_size(refseq_accession_number)
        partitions, regions = self._set_default_values(partitions, regions)

        data = self.rm_repeats_service.extract_info_by_chromosome(refseq_accession_number)

        self._graph_rm_options(graphs_config, data, size, partitions, regions, plot_type, save, dir, filename)


    @DBConnection
    @TryExcept
    @Timer
    def graph_rm_results_from_files_in_folder(self, directory_path: str, partitions: int, regions: int, plot_type: str, save: bool,
                                              dir: str):
        apply_function_to_files_in_folder(directory_path, graph_rm_results_from_file, partitions, regions,
                                          plot_type, save, dir)

    @DBConnection
    @TryExcept
    @Timer
    def graph_rm_results_of_genome_from_database(self, GCF: str, partitions: int, regions: int, plot_type: str, save: bool,
                                                 dir: str):
        chromosomes_ran_list = self.organisms_service.extract_chromosomes_refseq_accession_numbers_by_GCF(GCF)

        for refseq_accession_number in chromosomes_ran_list:
            self.graph_rm_results_from_database(refseq_accession_number, partitions, regions, plot_type, save, dir)

    @DBConnection
    @TryExcept
    @Timer
    def graph_genome_repeats_from_file(self, path: str, dir: str, partitions: int, regions: int, plot_type: str,
                                       save: bool):
        config = self.load_config()
        graphs_config = config.get('repeats', {})

        df = FileReader.read_repeats_results_file(path)
        df_list = FileReader.divide_genome_df_rows_by_chromosome(df)

        for index, df in enumerate(df_list):
            refseq_accession_number = df['refseq_accession_number'][0]
            filename, size = self.whole_chromosomes_service.extract_filename_and_size_by_refseq_accession_number(
                refseq_accession_number)

            partitions = int(partitions) if isinstance(partitions, str) else DEFAULT_PARTITIONS
            regions = int(regions) if isinstance(regions, str) else DEFAULT_REGIONS
            plot_type = plot_type or "line"

            logger.info(f"Starting generating graphs for chromosome {index + 1} - {refseq_accession_number}")

            # Call the refactored method to handle the graphing options
            self._graph_rm_options(graphs_config, df, size, partitions, regions, plot_type, save, dir, filename)

            logger.info(f"Completed the generation of graphs for chromosome {index + 1} - {refseq_accession_number}")

    @DBConnection
    @TryExcept
    @Timer
    def graph_recursive_from_database(self, refseq_accession_number: str, save: bool, name: str, n_max: int):
        config = self.load_config()
        graphs_config = config.get('recursive_repeats', {})

        n_max = n_max and int(n_max)
        filename = self.whole_chromosomes_service.extract_filename_by_refseq_accession_number(refseq_accession_number)
        data = self.recursive_repeats_whole_chromosomes_service.extract_info_by_chromosome(refseq_accession_number)

        if graphs_config.get('recursive_repeats_largest_values', False):
            Graphs.graph_recursive_repeats_largest_values_from_database(
                data, col="name", n_max=n_max, save=save, name=name, filename=filename)

        if graphs_config.get('grouped_by_recursive_repeat_length', False):
            Graphs.graph_grouped_by_recursive_repeat_length(
                data, col="name", save=save, name=name, filename=filename)

        if graphs_config.get('individual_plots_by_recursive_repeat_length', False):
            Graphs.graph_individual_plots_by_recursive_repeat_length(
                data, col="name", save=save, name=name, filename=filename)

    @DBConnection
    @TryExcept
    @Timer
    def graph_recursive_genome_from_database(self, GCF: str, save: bool, dir: str, n_max: int):
        n_max = n_max and int(n_max)
        refseq_accession_numbers = self.organism_service.extract_chromosomes_refseq_accession_numbers_by_GCF(GCF)

        for refseq_accession_number in refseq_accession_numbers:
            self.graph_recursive_from_database(refseq_accession_number, save, dir, n_max)

    @DBConnection
    @TryExcept
    @Timer
    def graph_gtf_from_file(self, path: str, partitions: int, regions: int, plot_type: str, save: bool, dir: str):
        config = self.load_config()
        graphs_config = config.get('genes', {})

        partitions = int(partitions) if isinstance(partitions, str) else DEFAULT_PARTITIONS
        regions = int(regions) if isinstance(regions, str) else DEFAULT_REGIONS
        plot_type = plot_type or "line"

        file_df = FileReader.read_gtf_file(path)
        chromosomes_df_list = FileReader.divide_genome_df_rows_by_chromosome(file_df)

        for df in chromosomes_df_list:
            refseq_accession_number = df['refseq_accession_number'][0]
            logger.info(f"Graphing for the sequence {refseq_accession_number}")
            try:
                chromosome_name, size = self.whole_chromosomes_service.extract_filename_and_size_by_refseq_accession_number(
                    refseq_accession_number)
            except Exception as e:
                logger.error(
                    f"Failed to extract whole chromosome ID for refseq_accession_number {refseq_accession_number}: {e}")
                break

            if graphs_config.get('distribution_of_genes_merged', False):
                Graphs.graph_distribution_of_genes_merged(
                    df, dir, size, partitions, regions, plot_type, chromosome_name, bool(save)
                )

            if graphs_config.get('distribution_of_genes', False):
                Graphs.graph_distribution_of_genes(
                    df, dir, legend=True, plot_type=plot_type, limit=20, regions=regions,
                    chromosome_name=chromosome_name, save=bool(save)
                )


    def _graph_gtf_single_chromosome_from_database(self, refseq_accession_number: str, name: str, partitions: int, regions: int,
                                                   plot_type: str, save: bool):
        config = self.load_config()
        graphs_config = config.get('genes', {})

        logger.info(f"Graphing for the sequence {refseq_accession_number}")
        chromosome_name, size = self.whole_chromosomes_service.extract_filename_and_size_by_refseq_accession_number(
                                    refseq_accession_number)
        df = self.gtf_genes_service.extract_genes_by_chromosome(refseq_accession_number)

        if graphs_config.get('distribution_of_genes_merged', False):
            Graphs.graph_distribution_of_genes_merged(df, name, size, partitions, regions, plot_type, chromosome_name,
                                                  bool(save))


    @DBConnection
    @TryExcept
    @Timer
    def graph_gtf_from_database(self, GCF: str, refseq_accession_number: str, partitions: int, regions: int, plot_type: str,
                                save: bool, dir: str):
        partitions = int(partitions) if isinstance(partitions, str) else DEFAULT_PARTITIONS
        regions = int(regions) if isinstance(regions, str) else DEFAULT_REGIONS
        plot_type = plot_type or "line"

        if GCF:
            if refseq_accession_number:
               self._graph_gtf_single_chromosome_from_database(refseq_accession_number, dir, partitions, regions, plot_type,
                                                          save)
            else:
                chromosomes_ran_list = self.organisms_service.extract_chromosomes_refseq_accession_numbers_by_GCF(GCF)
                logger.info(f"Graphing for the genome: {GCF}")
                for refseq_accession_number in chromosomes_ran_list:
                    self._graph_gtf_single_chromosome_from_database(refseq_accession_number, dir, partitions, regions,
                                                               plot_type, save)
        elif refseq_accession_number:
            self._graph_gtf_single_chromosome_from_database(refseq_accession_number, dir, partitions, regions, plot_type,
                                                       save)
        else:
            logger.error("Specify whether a GCF (organism/genome) or a refseq accession number (sequence/chromosome)")
            return



