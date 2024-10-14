import ast
from collections import Counter

import pandas as pd

from Loader import Loader
from src.Biocode.graphs.Graphs import Graphs
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.sequences.Sequence import Sequence
from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.GenesContainingRepeatsService import GenesContainingRepeatsService
from src.Biocode.services.LinearRepeatsWholeChromosomesService import LinearRepeatsWholeChromosomesService
from utils.FileReader import FileReader
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.folder import apply_function_to_files_in_folder
from utils.logger import logger


@Inject(repeats_service=RepeatsService,
        rm_repeats_service=RMRepeatsWholeChromosomesService,
        gtf_genes_service=GtfGenesService,
        whole_chromosomes_service=WholeChromosomesService,
        linear_repeats_whole_chromosomes_service=LinearRepeatsWholeChromosomesService,
        genes_containing_repeats_service=GenesContainingRepeatsService,
        file_reader=FileReader,
        loader=Loader)
class RepeatsLoader:

    def __init__(self, repeats_service=None, rm_repeats_service=None, whole_chromosomes_service=None, file_reader=None,
                 gtf_genes_service=None, linear_repeats_whole_chromosomes_service=None,
                 genes_containing_repeats_service=None, loader=None):
        self.repeats_service = repeats_service
        self.rm_repeats_service = rm_repeats_service
        self.gtf_genes_service = gtf_genes_service
        self.whole_chromosomes_service = whole_chromosomes_service
        self.linear_repeats_whole_chromosomes_service = linear_repeats_whole_chromosomes_service
        self.genes_containing_repeats_service = genes_containing_repeats_service
        self.file_reader = file_reader
        self.loader = loader

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
                                                    window_length=int(args.window_length))
            elif args.method == 'rm':
                logger.warning("Feature not implemented yet")

    def _find_kmers_recursively_in_genome(self, organism_name, gcf, data, k_range, save_to_db):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        genome_manager.find_only_kmers_recursively_and_calculate_multifractal_analysis_values(
            GCF=gcf, save_to_db=save_to_db, method_to_find_it="Recursively", k_range=k_range)

    def _find_kmers_linearly_in_genome(self, organism_name, gcf, data, k_range, save_to_db, dir, window_length):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        for manager in genome_manager.get_managers():
            self._find_kmers_linearly_in_sequence(gcf=gcf, sequence=manager.get_sequence(), k_range=k_range,
                                                  save_to_db=save_to_db, dir=dir, window_length=window_length)

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
                                                      dir=args.dir, window_length=int(args.window_length))
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
    def find_kmers_linearly_genes_sequence_command(self, args):
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

            self._find_kmers_linearly_genes_sequence(gcf=self.loader.get_gcf(), sequence=sequence,
                                                     k_range=ast.literal_eval(args.k_range), save_to_db=save_to_db,
                                                     dir=args.dir, window_length=int(args.window_length),
                                                     refseq_accession_number=self.loader.extract_refseq_accession_number(
                                                         args.path))
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _find_kmers_linearly_genes_sequence(self, gcf: str, sequence: Sequence, k_range, window_length,
                                            save_to_db: bool,
                                            dir: str, refseq_accession_number: str):

        sequence_manager = RegionSequenceManager(sequence=sequence)
        sequence_nts = sequence_manager.get_sequence().get_sequence()

        genes_sequence_df = self.gtf_genes_service.extract_genes_by_chromosome(refseq_accession_number)

        # CHECK: TRY TO PERFORM THIS WITHOUT CALCULATING window_profiles
        window_profiles, most_frequent_nplets = self.__get_kmers_mapped_in_windows(sequence_manager, k_range,
                                                                                   window_length)
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

        self.load_linear_repeats(window_profiles_only_genes, refseq_accession_number, most_frequent_nplets)
        self.load_genes_containing_repeats(refseq_accession_number, sequence_nts)

        # Graphs.plot_combined_kmer_frequency(window_profiles_only_genes, most_frequent_nplets, sequence_manager.get_sequence_name(),
        #                                        dir, True, subfolder="linear_repeats_genes")

        Graphs.plot_combined_kmer_frequency_graph_per_k(window_profiles_only_genes, most_frequent_nplets,
                                                        sequence_manager.get_sequence_name(), dir, True,
                                                        subfolder="linear_repeats_genes")

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

            # Count overlapping occurrences of the repeat in the gene sequence - For non overlapping use count method
            total_count += self.___count_overlapping_occurrences(gene_sequence, repeat)

        return total_count

    def ___count_overlapping_occurrences(self, sequence: str, repeat: str) -> int:
        """
        Count overlapping occurrences of a substring (repeat) in a given sequence.

        Parameters:
        - sequence (str): The string to search within.
        - repeat (str): The substring to count occurrences of.

        Returns:
        - int: The number of overlapping occurrences of the repeat in the sequence.
        """
        count = start = 0
        while True:
            start = sequence.find(repeat, start)
            if start == -1:  # No more occurrences found
                break
            count += 1
            start += 1  # Move one character forward to allow for overlaps
        return count


    def _find_kmers_linearly_in_sequence(self, gcf: str, sequence: Sequence, k_range: tuple, save_to_db: bool, dir: str,
                                         window_length: int):
        sequence_manager = RegionSequenceManager(sequence=sequence)

        #self._plot_all_kmers(window_profiles, most_frequent_nplets, sequence_manager.get_sequence_name(), save_to_db, dir)

        window_profiles, most_frequent_nplets = self.__get_kmers_mapped_in_windows(sequence_manager, k_range, window_length)
        Graphs.plot_combined_kmer_frequency(window_profiles, most_frequent_nplets, sequence_manager.get_sequence_name(),
                                            dir, True, subfolder="linear_repeats_all")
        Graphs.plot_combined_kmer_frequency_graph_per_k(window_profiles, most_frequent_nplets, sequence_manager.get_sequence_name(),
                                                        dir, True, subfolder="linear_repeats_all/per_k")
        if save_to_db:
            logger.warn("window_profiles")
            logger.warn(window_profiles)
            logger.warn("******************")
            logger.warn(most_frequent_nplets)
            self.load_linear_repeats(window_profiles, refseq_accession_number=sequence.get_refseq_accession_number(),
                                     most_frequent_nplets=most_frequent_nplets)


    def __get_kmers_mapped_in_windows(self, sequence_manager: RegionSequenceManager, k_range: tuple, window_length: int):
        sequence = sequence_manager.get_sequence().get_sequence()
        result = {}
        for k in range(k_range[0], k_range[1] + 1):
            result[f'{k}-mers'] = self.__count_kmers(sequence, k)

        most_frequent_nplets = self.__get_most_frequent_nplets(nplet_counts=result, top_n=10)
        window_profiles = self.__map_kmers_in_windows(sequence, most_frequent_nplets, window_length=window_length)
        return window_profiles, most_frequent_nplets


    def __get_most_frequent_nplets(self, nplet_counts, top_n=10):
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
                window_count[k] = {kmer: self.__count_kmers_in_sequence(window_seq, kmer) for kmer, _ in kmers}
            window_profiles.append(window_count)
        return window_profiles

    def __count_kmers_in_sequence(self, sequence, kmer):
        count = 0
        k = len(kmer)

        # Slide through the sequence and count overlapping occurrences of `kmer`
        for i in range(len(sequence) - k + 1):
            if sequence[i:i + k] == kmer:
                count += 1
        return count

    def __count_kmers(self, sequence, k):
        kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
        return Counter(kmers)

    def _plot_all_kmers(self, window_profiles, most_frequent_nplets, sequence_name, save_to_db, dir):
        """
        Plots the frequency profile for all most frequent k-mers in the genome.
        """
        for k, kmers in most_frequent_nplets.items():
            for kmer, _ in kmers:
                Graphs.plot_kmer_frequency(window_profiles, kmer, sequence_name, dir, True)

    def load_linear_repeats(self, window_profiles: list[dict], refseq_accession_number: str,
                            most_frequent_nplets: dict):
        for window in window_profiles:
            for kmer_length, kmer_dict in window.items():
                for kmer, count in kmer_dict.items():
                    record = (kmer, '', 'Linear')
                    repeat_id = self.repeats_service.insert(record=record)
                    whole_chromosome_id = self.whole_chromosomes_service.extract_id_by_refseq_accession_number(
                        refseq_accession_number)

                    # Check if kmer_length exists in most_frequent_nplets
                    if kmer_length in most_frequent_nplets:
                        nplet_list = most_frequent_nplets[kmer_length]
                        nplet_dict = dict(nplet_list)
                        repeat_count = nplet_dict.get(kmer, 0) # zero if not found
                    else:
                        repeat_count = 0

                    self.linear_repeats_whole_chromosomes_service.insert(record=(
                        repeat_id, whole_chromosome_id, len(kmer), repeat_count
                    ))

    def load_genes_containing_repeats(self, refseq_accession_number: str, sequence_nts: str):
        genes_df = self.gtf_genes_service.extract_genes_by_chromosome(refseq_accession_number)
        repeats_df = self.linear_repeats_whole_chromosomes_service.extract_linear_repeats_by_refseq_accession_number(
            refseq_accession_number)

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

                # Use the existing method to count overlapping occurrences
                count = self.___count_overlapping_occurrences(gene_sequence, repeat_sequence)

                # If count > 0, save the record
                if count > 0:
                    record = (gene_id, repeat_id, count)
                    self.genes_containing_repeats_service.insert(record=record)

    def _extract_gene_sequence_based_on_positions(self, sequence_nts: str, start_position: int, end_position: int):
        return sequence_nts[start_position:end_position + 1]
