import ast
from collections import Counter

import pandas as pd

from Loader import Loader
from src.Biocode.graphs.Graphs import Graphs
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.sequences.Sequence import Sequence
from src.Biocode.services.GenesContainingRepeatsService import GenesContainingRepeatsService
from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.services.LinearRepeatsRegionChromosomesService import LinearRepeatsRegionChromosomesService
from src.Biocode.services.LinearRepeatsWholeChromosomesService import LinearRepeatsWholeChromosomesService
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RegionChromosomesService import RegionChromosomesService
from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.utils.utils import (adapt_dataframe_to_window_profiles, adapt_dataframe_to_most_frequent_nplets)
from utils.FileReader import FileReader
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.folder import apply_function_to_files_in_folder
from utils.logger import logger


@Inject(repeats_service=RepeatsService,
        rm_repeats_service=RMRepeatsWholeChromosomesService,
        gtf_genes_service=GtfGenesService,
        organisms_service=OrganismsService,
        whole_chromosomes_service=WholeChromosomesService,
        region_chromosomes_service=RegionChromosomesService,
        linear_repeats_whole_chromosomes_service=LinearRepeatsWholeChromosomesService,
        linear_repeats_region_chromosomes_service=LinearRepeatsRegionChromosomesService,
        genes_containing_repeats_service=GenesContainingRepeatsService,
        file_reader=FileReader,
        loader=Loader)
class RepeatsLoader:

    def __init__(self, repeats_service=None, rm_repeats_service=None, whole_chromosomes_service=None,
                 region_chromosomes_service=None, file_reader=None, organisms_service=None,
                 gtf_genes_service=None, linear_repeats_whole_chromosomes_service=None,
                 linear_repeats_region_chromosomes_service=None,
                 genes_containing_repeats_service=None, loader=None):
        self.repeats_service = repeats_service
        self.rm_repeats_service = rm_repeats_service
        self.gtf_genes_service = gtf_genes_service
        self.organisms_service = organisms_service
        self.whole_chromosomes_service = whole_chromosomes_service
        self.region_chromosomes_service = region_chromosomes_service
        self.linear_repeats_whole_chromosomes_service = linear_repeats_whole_chromosomes_service
        self.linear_repeats_region_chromosomes_service = linear_repeats_region_chromosomes_service
        self.genes_containing_repeats_service = genes_containing_repeats_service
        self.file_reader = file_reader
        self.loader = loader
        self.IS_COVERED = False

    @DBConnection
    @TryExcept
    @Timer
    def load_RM_repeats_from_file(self, path):
        logger.info(f"Loading repeats results from RepeatMasker from {path}")
        df = self.file_reader.read_repeats_results_file(path)
        self._load_repeats_file_to_database(df, 'RM')

    @DBConnection
    @TryExcept
    @Timer
    def load_RM_repeats_from_folder(self, path):
        apply_function_to_files_in_folder(path, self.load_RM_repeats_from_file)

    @DBConnection
    @TryExcept
    @Timer
    def load_genome_repeats_file(self, path):
        logger.info(f"Loading repeats results of a genome from a single file in {path}")
        df = self.file_reader.read_repeats_results_file(path)
        df_list = self.file_reader.divide_genome_df_rows_by_chromosome(df)

        for df in df_list:
            self._load_repeats_file_to_database(df, 'Genome in file')

    def _load_repeats_file_to_database(self, df: pd.DataFrame, method_to_find_it: str):
        logger.info(f"*********** Starting loading the chromosome: {df['refseq_accession_number'][0]} ***********")

        for _, row in df.iterrows():
            repeat_id = self.repeats_service.insert(record=(row['repeat'], row['class_family'], method_to_find_it))
            whole_chromosome_id = self._get_whole_chromosome_id(row['refseq_accession_number'])
            if whole_chromosome_id:
                self._insert_rm_repeat(row, repeat_id, whole_chromosome_id)

        logger.info(f"*********** Completed loading the chromosome: {df['refseq_accession_number'][0]} ***********")

    def _get_whole_chromosome_id(self, refseq_accession_number):
        try:
            return self.whole_chromosomes_service.extract_id_by_refseq_accession_number(refseq_accession_number)
        except Exception as e:
            logger.error(
                f"Failed to extract whole chromosome ID for refseq_accession_number {refseq_accession_number}: {e}")
            return None

    def _insert_rm_repeat(self, row, repeat_id, whole_chromosome_id):
        record = (
            repeat_id,
            whole_chromosome_id,
            row['sw_score'],
            row['percentage_divergence'],
            row['percentage_deletions'],
            row['percentage_insertions'],
            row['query_begin'],
            row['query_end'],
            row['repeat_length'],
            row['query_left'],
            row['strand'],
            row['repeat_begin'],
            row['repeat_end'],
            row['repeat_left']
        )
        self.rm_repeats_service.insert(record=record)

    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_genome_command(self, args):
        save_to_db = False if args.save_to_db == 'false' else True
        self.organism = args.name
        self.loader.set_organism(self.organism)
        #self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
        #                   amount_chromosomes=self.loader.get_amount_chromosomes())
        if args.method:
            if args.method == 'r':

                self._find_kmers_recursively_in_genome(organism_name=self.loader.get_organism_name(),
                                                       gcf=self.loader.get_gcf(),
                                                       data=self.loader.get_data(),
                                                       k_range=ast.literal_eval(args.k_range),
                                                       save_to_db=save_to_db)
            elif args.method == "l":
                self._find_kmers_linearly_in_genome(organism_name=self.loader.get_organism_name(),
                                                    gcf=self.loader.get_gcf(),
                                                    data=self.loader.get_data(), k_range=ast.literal_eval(args.k_range),
                                                    save_to_db=save_to_db, dir=args.dir,
                                                    window_length=int(args.window_length),
                                                    graph_from_file=bool(args.graph_from_file))
            elif args.method == 'rm':
                logger.warning("Feature not implemented yet")

    def _find_kmers_recursively_in_genome(self, organism_name, gcf, data, k_range, save_to_db):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        genome_manager.find_only_kmers_recursively_and_calculate_multifractal_analysis_values(
            GCF=gcf, save_to_db=save_to_db, method_to_find_it="Recursively", k_range=k_range)

    def _find_kmers_linearly_in_genome(self, organism_name, gcf, data, k_range, save_to_db, dir, window_length,
                                       graph_from_file):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        for manager in genome_manager.get_managers():
            self._find_kmers_linearly_in_sequence(gcf=gcf, sequence=manager.get_sequence(), k_range=k_range,
                                                  save_to_db=save_to_db, dir=dir, window_length=window_length,
                                                  graph_from_file=graph_from_file)

    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_sequence_command(self, args):
        if args.name:
            self.organism = args.name
            self.loader.set_organism(self.organism)
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

        if args.path:
            sequence = Sequence(sequence=self.loader.read_fasta_sequence(file_path=args.path),
                                name=self.loader.extract_file_name(file_path=args.path),
                                organism_name=self.loader.get_organism_name(),
                                refseq_accession_number=self.loader.extract_refseq_accession_number(args.path))
            save_to_db = False if args.save_to_db == 'false' else True

            #self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
            #                    amount_chromosomes=self.loader.get_amount_chromosomes())

            if args.method == 'r':
                self._find_kmers_recursively_in_sequence(gcf=self.loader.get_gcf(), sequence=sequence,
                                                         k_range=ast.literal_eval(args.k_range), save_to_db=save_to_db)
            elif args.method == 'l':
                self._find_kmers_linearly_in_sequence(gcf=self.loader.get_gcf(), sequence=sequence,
                                                      k_range=ast.literal_eval(args.k_range), save_to_db=save_to_db,
                                                      dir=args.dir, window_length=int(args.window_length),
                                                      graph_from_file=bool(args.graph_from_file))
            elif args.method == 'rm':
                logger.warning("Feature not implemented yet")
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _find_kmers_recursively_in_sequence(self, gcf, sequence, k_range, save_to_db):
        sequence_manager = SequenceManager(sequence=sequence)
        sequence_manager.calculate_multifractal_analysis_values(gcf)
        kmers_list = sequence_manager.find_only_kmers_recursively(k_range)
        if save_to_db:
            sequence_manager.save_repeats_found_recursively_to_db(
                kmers_list=kmers_list, GCF=gcf, method_to_find_it="Recursively")

    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_linearly_genes_genome_command(self, args):
        if args.name:
            self.organism = args.name
            self.loader.set_organism(self.organism)
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")
        apply_function_to_files_in_folder(self.loader.get_organism_path(),
                                          self.find_kmers_linearly_genes_sequence_command, args)


    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_linearly_genes_sequence_command(self, filepath: str, args):
        if args.name:
            self.organism = args.name
            self.loader.set_organism(self.organism)
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

        if filepath:
            sequence = Sequence(sequence=self.loader.read_fasta_sequence(file_path=filepath),
                                name=self.loader.extract_file_name(file_path=filepath),
                                organism_name=self.loader.get_organism_name(),
                                refseq_accession_number=self.loader.extract_refseq_accession_number(filepath))
            save_to_db = False if args.save_to_db == 'false' else True

            #self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
            #                    amount_chromosomes=self.loader.get_amount_chromosomes())

            self._find_kmers_linearly_genes_sequence(sequence=sequence, save_to_db=save_to_db, dir=args.dir,
                                                     window_length=int(args.window_length),
                                                     refseq_accession_number=self.loader.extract_refseq_accession_number(
                                                         filepath), graph_from_file=bool(args.graph_from_file),
                                                     size=int(args.size))
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _find_kmers_linearly_genes_sequence(self, sequence: Sequence, window_length,
                                            save_to_db: bool, dir: str, refseq_accession_number: str,
                                            graph_from_file: bool, size: int):
        sequence_manager = RegionSequenceManager(sequence=sequence, window_length=window_length)
        sequence_nts = sequence_manager.get_sequence().get_sequence()

        genes_sequence_df = self.gtf_genes_service.extract_genes_by_chromosome(refseq_accession_number)
        logger.info("Gtf_genes dataframe of chromosome extracted from database")
        #window_profiles, most_frequent_nplets = self.__get_kmers_mapped_in_windows(sequence_manager,k_range,window_length)
        window_profiles, most_frequent_nplets = self.__get_kmers_mapped_in_windows_from_database(
            refseq_accession_number)
        logger.info(
            "RepeatsLoader._find_kmers_linearly_genes_sequence: window_profiles and most_frequent_nplets extracted from database")

        logger.info("Starting building window_profiles_only_genes")
        window_profiles_only_genes = []
        for window_index, window in enumerate(window_profiles):
            kmer_key_dict = {}
            for kmer_key, repeats_dict_value in window.items():  # (4-12) mers
                repeat_count_dict = {}
                for repeat_key, count_value in repeats_dict_value.items():  # 10 repeats
                    repeat_count_in_genes = self.__count_repeat_in_genes_given_a_sequence_chunk(
                        genes_sequence_df,
                        repeat_key,
                        sequence_nts,
                        window_index * window_length,
                        (window_index + 1) * window_length - 1)
                    repeat_count_dict[repeat_key] = repeat_count_in_genes
                kmer_key_dict[kmer_key] = repeat_count_dict
            window_profiles_only_genes.append(kmer_key_dict)
        logger.info("Completed building of window_profile_only_genes")
        if save_to_db:
            self.load_linear_repeats(window_profiles_only_genes, refseq_accession_number, most_frequent_nplets,
                                     sequence_manager.get_regions_refseq_accessions_numbers(),
                                     method_to_find_it='Linear in genes')
            logger.info("Completed loading of 'Linear in genes' repeats")
            self.load_genes_containing_repeats(refseq_accession_number, sequence_nts, genes_sequence_df, size)
            logger.info("Completed loading of genes_containing_repeats")

        # Graphs.plot_combined_kmer_frequency(window_profiles_only_genes, most_frequent_nplets, sequence_manager.get_sequence_name(),
        #                                        dir, True, window_length, subfolder="linear_repeats_genes")

        if graph_from_file:
            Graphs.plot_combined_kmer_frequency_graph_per_k(window_profiles_only_genes, most_frequent_nplets,
                                                            sequence_manager.get_sequence_name(), dir, True,
                                                            window_length,
                                                            subfolder=f"linear_repeats_genes/per_k/{sequence.get_name()}")

    def __get_kmers_mapped_in_windows_from_database(self, refseq_accession_number: str):
        whole_repeats_df = self.linear_repeats_whole_chromosomes_service.extract_linear_repeats_by_refseq_accession_number(
            refseq_accession_number)
        region_repeats_df = self.linear_repeats_region_chromosomes_service.extract_linear_repeats_by_refseq_accession_number(
            refseq_accession_number)

        window_length = region_repeats_df.iloc[0]['window_length']
        window_profiles = adapt_dataframe_to_window_profiles(region_repeats_df)
        most_frequent_nplets = adapt_dataframe_to_most_frequent_nplets(whole_repeats_df)
        return window_profiles, most_frequent_nplets

    def __count_repeat_in_genes_given_a_sequence_chunk(self, genes_df, repeat: str, sequence: str,
                                                       initial_position_value, final_position_value):
        """
        Count the occurrences of a specific repeat in the sequences of genes
        that fall within a specified range of positions.

        Parameters:
        - genes_df (pd.DataFrame): DataFrame containing gene information with columns "start_position", "end_position", and "strand".
        - repeat (str): The nucleotide repeat to search for.
        - sequence (str): The full nucleotide sequence to analyze.
        - initial_position_value (int): The starting position of the chunk to analyze (1-based).
        - final_position_value (int): The ending position of the chunk to analyze (1-based).

        Returns:
        - int: The total number of appearances of the repeat in the specified gene segments.
        """

        total_count = 0

        for index, row in genes_df.iterrows():
            start = row['start_position'] - 1  # Convert to 0-based index
            end = row['end_position'] - 1  # Convert to 0-based index

            # Check if the gene is within the chunk
            if end < initial_position_value:
                continue  # Gene ends before the chunk starts
            if start > final_position_value:
                break  # Gene starts after the chunk ends

            # Calculate the overlapping range for this gene within the chunk
            overlap_start = max(start, initial_position_value)
            overlap_end = min(end, final_position_value)

            # Extract the relevant segment of the sequence
            gene_sequence = sequence[overlap_start:overlap_end + 1]

            total_count += self.__count_non_overlapping_kmers_having_kmer(gene_sequence, repeat)

        return total_count

    def _find_kmers_linearly_in_sequence(self, gcf: str, sequence: Sequence, k_range: tuple, save_to_db: bool, dir: str,
                                         window_length: int, graph_from_file: bool):
        sequence_manager = RegionSequenceManager(sequence=sequence, window_length=window_length)
        logger.info(
            f"Starting analysis for finding kmers and window profiles of sequence {sequence.get_name()} - {sequence.get_organism_name()}")
        window_profiles, most_frequent_nplets = self.__get_kmers_mapped_in_windows(sequence_manager, k_range,
                                                                                   window_length)
        logger.info(
            f"Analysis completed for finding kmers and window profiles of sequence {sequence.get_name()} - {sequence.get_organism_name()}")

        if graph_from_file:
            # self._plot_all_kmers(window_profiles, most_frequent_nplets, sequence_manager.get_sequence_name(), save_to_db, dir)

            Graphs.plot_combined_kmer_frequency(window_profiles, most_frequent_nplets,
                                                sequence_manager.get_sequence_name(),
                                                dir, True, window_length, subfolder="linear_repeats_all")
            Graphs.plot_combined_kmer_frequency_graph_per_k(window_profiles, most_frequent_nplets,
                                                            sequence_manager.get_sequence_name(),
                                                            dir, True, window_length,
                                                            subfolder=f"linear_repeats_all/per_k/{sequence.get_name()}")
            logger.info(f"Graphing plot from file completed for {sequence.get_name()} - {sequence.get_organism_name()}")
        if save_to_db:
            logger.info(f"Loading repeats/kmers found linearly to database for {sequence.get_name()}")
            self.load_linear_repeats(window_profiles, refseq_accession_number=sequence.get_refseq_accession_number(),
                                     most_frequent_nplets=most_frequent_nplets,
                                     regions_refseq_accession_number_list=sequence_manager.get_regions_refseq_accessions_numbers(),
                                     method_to_find_it='Linear')

    def __get_kmers_mapped_in_windows(self, sequence_manager: RegionSequenceManager, k_range: tuple,
                                      window_length: int):
        sequence = sequence_manager.get_sequence().get_sequence()
        result = {}
        for k in range(k_range[0], k_range[1] + 1):
            result[f'{k}-mers'] = self.__count_overlapping_kmers_having_k(sequence, k)

        most_frequent_nplets = self.__get_most_frequent_nplets(nplet_counts=result, top_n=10)
        window_profiles = self.__map_kmers_in_windows(sequence, most_frequent_nplets, window_length=window_length)
        return window_profiles, most_frequent_nplets

    def __get_most_frequent_nplets(self, nplet_counts, top_n=10) -> dict[str, list[tuple[str, int]]]:
        most_frequent = {}
        for k, counts in nplet_counts.items():
            most_frequent[k] = counts.most_common(top_n)
        return most_frequent

    def __map_kmers_in_windows(self, genome_sequence: str, k_mers, window_length):
        window_profiles = []

        for start in range(0, len(genome_sequence), window_length):
            window_seq = genome_sequence[start:start + window_length]
            window_count = {}
            for k, kmers in k_mers.items():
                window_count[k] = {kmer: self.__count_non_overlapping_kmers_having_kmer(window_seq, kmer) for kmer, _ in
                                   kmers}
            window_profiles.append(window_count)
        return window_profiles

    def __count_non_overlapping_kmers_having_kmer(self, sequence: str, kmer: str) -> int:
        def __aux_count_non_overlapping(seq: str, kmer: str) -> int:
            count = 0
            start = 0

            while True:
                start = seq.find(kmer, start)
                if start == -1:
                    break
                count += 1
                start += len(kmer)

            return count

        if self.IS_COVERED:
            return __aux_count_non_overlapping(sequence, kmer)
        else:
            valid_segments = self.___divide_in_valid_segments_of_nucleotides(sequence)
            return sum(__aux_count_non_overlapping(segment, kmer) for segment in valid_segments)

    def __count_overlapping_kmers_having_kmer(self, sequence: str, kmer: str) -> int:
        k = len(kmer)
        count = 0

        if self.IS_COVERED:
            for i in range(len(sequence) - k + 1):
                if sequence[i:i + k] == kmer:
                    count += 1
        else:
            valid_segments = self.___divide_in_valid_segments_of_nucleotides(sequence)
            for segment in valid_segments:
                count += sum(1 for i in range(len(segment) - k + 1) if segment[i:i + k] == kmer)

        return count

    def __count_overlapping_kmers_with_non_overlapping_counts_having_k(self, sequence, k):
        overlapping_counts = self.__count_overlapping_kmers_having_k(sequence, k)
        non_overlapping_counts = Counter()

        # Count non-overlapping occurrences for each k-mer
        for kmer in overlapping_counts.keys():
            count = self.__count_non_overlapping_kmers_having_kmer(sequence, kmer)
            non_overlapping_counts[kmer] = count

        return non_overlapping_counts

    def __count_non_overlapping_kmers_having_k(self, sequence: str, k: int) -> Counter:
        if self.IS_COVERED:
            kmers = [sequence[i:i + k] for i in range(0, len(sequence) - k + 1, k)]
        else:
            valid_segments = self.___divide_in_valid_segments_of_nucleotides(sequence)
            kmers = [
                segment[i:i + k] for segment in valid_segments
                for i in range(0, len(segment) - k + 1, k)
            ]

        return Counter(kmers)

    def ___divide_in_valid_segments_of_nucleotides(self, sequence: str) -> list[str]:
        valid_segments = []
        current_segment = []

        for nuc in sequence:
            if nuc in 'ATCG':
                current_segment.append(nuc)
            else:
                if current_segment:
                    valid_segments.append(''.join(current_segment))
                    current_segment = []

        # Add the last segment if it exists
        if current_segment:
            valid_segments.append(''.join(current_segment))
        return valid_segments

    def __count_overlapping_kmers_having_k(self, sequence: str, k: int) -> Counter:
        if self.IS_COVERED:
            kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
        else:
            valid_segments = self.___divide_in_valid_segments_of_nucleotides(sequence)
            kmers = [segment[i:i + k] for segment in valid_segments for i in range(len(segment) - k + 1)]

        return Counter(kmers)

    def _plot_all_kmers(self, window_profiles, most_frequent_nplets, sequence_name, save, dir):
        """
        Plots the frequency profile for all most frequent k-mers in the genome.
        """
        for k, kmers in most_frequent_nplets.items():
            for kmer, _ in kmers:
                Graphs.plot_kmer_frequency(window_profiles, kmer, sequence_name, dir, save)

    def load_linear_repeats(self, window_profiles: list[dict], refseq_accession_number: str,
                            most_frequent_nplets: dict, regions_refseq_accession_number_list: list[str],
                            method_to_find_it: str):
        for window_index, window in enumerate(window_profiles):
            for kmer_length, kmer_dict in window.items():
                for kmer, count in kmer_dict.items():
                    record = (kmer, '', method_to_find_it)
                    repeat_id = self.repeats_service.insert(record=record)
                    whole_chromosome_id = self.whole_chromosomes_service.extract_id_by_refseq_accession_number(
                        refseq_accession_number)
                    region_chromosome_id = self.region_chromosomes_service.extract_id_by_refseq_accession_number(
                        regions_refseq_accession_number_list[window_index]
                    )

                    # Check if kmer_length exists in most_frequent_nplets
                    if kmer_length in most_frequent_nplets:
                        nplet_list = most_frequent_nplets[kmer_length]
                        nplet_dict = dict(nplet_list)
                        repeat_count_in_whole_chr = nplet_dict.get(kmer, 0)  # zero if not found
                    else:
                        repeat_count_in_whole_chr = 0

                    self.linear_repeats_whole_chromosomes_service.insert(record=(
                        repeat_id, whole_chromosome_id, len(kmer), repeat_count_in_whole_chr
                    ))
                    self.linear_repeats_region_chromosomes_service.insert(record=(
                        repeat_id, region_chromosome_id, len(kmer), count
                    ))

    def load_genes_containing_repeats(self, refseq_accession_number: str, sequence_nts: str, genes_df: pd.DataFrame,
                                      size: int):
        repeats_df = self.linear_repeats_whole_chromosomes_service.extract_linear_in_genes_repeats_by_refseq_accession_number(
            refseq_accession_number, k_range=(size, size))
        # Iterate over each gene in genes_df
        for index, gene in genes_df.iterrows():
            gene_id = gene['id']
            start_position = gene['start_position']
            end_position = gene['end_position']

            gene_sequence = self._extract_gene_sequence_based_on_positions(sequence_nts, start_position, end_position)

            # Iterate over each repeat in repeats_df
            for _, repeat in repeats_df.iterrows():
                repeat_id = repeat['repeats_id']
                repeat_sequence = repeat['name']

                count = self.__count_non_overlapping_kmers_having_kmer(gene_sequence, repeat_sequence)

                # If count > 0, save the record
                if count > 0:
                    record = (gene_id, repeat_id, count)
                    self.genes_containing_repeats_service.insert(record=record)

    def _extract_gene_sequence_based_on_positions(self, sequence_nts: str, start_position: int, end_position: int):
        return sequence_nts[start_position:end_position + 1]
