from Loader import Loader
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.sequences.Sequence import Sequence
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.logger import logger
from src.Biocode.graphs.Graphs import Graphs
from collections import Counter
import ast

@Inject(organisms_service = OrganismsService, loader = Loader, gtf_genes_service = GtfGenesService)
class Analyzer:
    def __init__(self, organisms_service: OrganismsService = None, loader: Loader = None,
                 gtf_genes_service: GtfGenesService = None):
        self.organisms_service = organisms_service
        self.gtf_genes_service = gtf_genes_service
        self.loader = loader
        self.organism = ""

    @Timer
    def _load_organism(self, organism_name, gcf, amount_chromosomes):
        logger.info(f"Loading organism {organism_name} - {gcf} - {amount_chromosomes} chromosomes - to database")
        self.organisms_service.insert(record=(organism_name, gcf, amount_chromosomes))

    @DBConnection
    @TryExcept
    @Timer
    def analyze_genome_command(self, args):
        save_to_db = False if args.save_to_db == 'false' else True

        if args.name:
            self.organism= args.name
            self.loader.set_organism(self.organism)
            self._validate_mode_analyzing_genome(args, save_to_db=save_to_db)

        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

    def _validate_mode_analyzing_genome(self, args, save_to_db: bool):
        if args.mode:
            if args.mode == 'whole':
                self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                              amount_chromosomes=self.loader.get_amount_chromosomes())
                self.__whole_MFA_genome(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(), data=self.loader.get_data(),
                                        save_to_db=save_to_db)
            elif args.mode == 'regions':
                self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                              amount_chromosomes=self.loader.get_amount_chromosomes())
                self.__regions_MFA_genome(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(), data=self.loader.get_data(),
                                          regions_number=args.regions_number, window_length=args.window_length,
                                          save_to_db=save_to_db)
            else:
                raise Exception("Enter a valid mode (whole or regions)")
        else:
            raise Exception("Enter a valid mode (whole or regions)")

    def __whole_MFA_genome(self, organism_name, gcf, data, save_to_db):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db)
        #genome_manager.graph_linear_fit()
        #genome_manager.generate_df_results()

    def __regions_MFA_genome(self, organism_name, gcf, data, regions_number, window_length, save_to_db):
        window_length = int(window_length) if window_length else None
        regions_number = int(regions_number) if regions_number else None
        region_genome_manager = RegionGenomeManager(genome_data=data, regions_number=regions_number,
                                                    window_length=window_length, organism_name=organism_name)
        region_genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db)

    @DBConnection
    @TryExcept
    @Timer
    def analyze_sequence_command(self, args):
        if args.name:
            self.organism= args.name
            self.loader.set_organism(self.organism)
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

        if args.path:
            sequence = Sequence(sequence=self.loader.read_fasta_sequence(file_path=args.path),
                                name=self.loader.extract_file_name(file_path=args.path),
                                organism_name=self.loader.get_organism_name(),
                                refseq_accession_number=self.loader.extract_refseq_accession_number(args.path))
            save_to_db = False if args.save_to_db == 'false' else True
            self._validate_mode_analyzing_sequence(args, sequence, save_to_db=save_to_db)
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _validate_mode_analyzing_sequence(self, args, sequence: Sequence, save_to_db: bool):
        if args.mode:
            if args.mode == 'whole':
                self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                              amount_chromosomes=self.loader.get_amount_chromosomes())
                self.__whole_MFA_sequence(gcf=self.loader.get_gcf(), sequence=sequence, save_to_db=save_to_db)
            elif args.mode == 'regions':
                self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                              amount_chromosomes=self.loader.get_amount_chromosomes())
                self.__regions_MFA_sequence(gcf=self.loader.get_gcf(), sequence=sequence,
                                     regions_number=args.regions_number, window_length=args.window_length,
                                            save_to_db=save_to_db)
            else:
                raise Exception("Enter a valid mode (whole or regions)")
        else:
            raise Exception("Enter a valid mode (whole or regions)")

    def __whole_MFA_sequence(self, gcf, sequence, save_to_db):
        sequence_manager = SequenceManager(sequence=sequence)
        sequence_manager.calculate_multifractal_analysis_values(gcf)

        if save_to_db:
            sequence_manager.save_results_to_db_during_execution(GCF=gcf)

    def __regions_MFA_sequence(self, gcf, sequence, regions_number, window_length, save_to_db):
        window_length = int(window_length) if window_length else None
        regions_number = int(regions_number) if regions_number else None
        # regions_number is None here
        region_sequence_manager = RegionSequenceManager(sequence=sequence, regions_number=regions_number,
                                                        window_length=window_length)
        region_sequence_manager.calculate_multifractal_analysis_values(gcf)

        if save_to_db:
            region_sequence_manager.save_results_to_db_during_execution(GCF=gcf)

    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_genome_command(self, args):
        save_to_db = False if args.save_to_db == 'false' else True
        self.organism = args.name
        self.loader.set_organism(self.organism)
        self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                            amount_chromosomes=self.loader.get_amount_chromosomes())
        if args.method:
            if args.method == 'r':

                self._find_kmers_recursively_in_genome(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                                                 data=self.loader.get_data(), k_range=ast.literal_eval(args.k_range),
                                                       save_to_db=save_to_db)
            elif args.method == "l":
                self._find_kmers_linearly_in_genome(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                                                    data=self.loader.get_data(), k_range=ast.literal_eval(args.k_range),
                                                    save_to_db=save_to_db, dir=args.dir, window_length=int(args.window_length))
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

            self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                                amount_chromosomes=self.loader.get_amount_chromosomes())

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


    def _find_kmers_linearly_in_sequence(self, gcf: str, sequence: Sequence, k_range: tuple, save_to_db: bool, dir: str,
                                         window_length: int):
        sequence_manager = RegionSequenceManager(sequence=sequence)

        #self._plot_all_kmers(window_profiles, most_frequent_nplets, sequence_manager.get_sequence_name(), save_to_db, dir)

        window_profiles, most_frequent_nplets = self.__get_kmers_mapped_in_windows(sequence_manager, k_range, window_length)
        Graphs.plot_combined_kmer_frequency(window_profiles, most_frequent_nplets, sequence_manager.get_sequence_name(),
                                            dir, True, subfolder="linear_repeats_all")

        if save_to_db:
            #sequence_manager.save_repeats_found_linearly_to_db(
            #)
            pass
            # MISSING HERE


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

            self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                                amount_chromosomes=self.loader.get_amount_chromosomes())

            self._find_kmers_linearly_genes_sequence(gcf=self.loader.get_gcf(), sequence=sequence,
                                                      k_range=ast.literal_eval(args.k_range), save_to_db=save_to_db,
                                                      dir=args.dir, window_length=int(args.window_length),
                                                     refseq_accession_number=self.loader.extract_refseq_accession_number(args.path))
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _find_kmers_linearly_genes_sequence(self, gcf: str, sequence: Sequence, k_range, window_length, save_to_db: bool,
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
            for kmer_key, repeats_dict_value in window.items(): # (4-12) mers
                repeat_count_dict = {}
                for repeat_key, count_value in repeats_dict_value.items(): # 10 repeats
                    repeat_count_in_genes = self.__count_repeat_in_genes_given_a_sequence_chunk(
                        genes_sequence_df,
                        repeat_key,
                        sequence_nts,
                        window_index * window_length,
                        (window_index + 1) * window_length - 1)
                    repeat_count_dict[repeat_key] = repeat_count_in_genes
                kmer_key_dict[kmer_key] = repeat_count_dict
            window_profiles_only_genes.append(kmer_key_dict)
        #Graphs.plot_combined_kmer_frequency(window_profiles_only_genes, most_frequent_nplets, sequence_manager.get_sequence_name(),
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


