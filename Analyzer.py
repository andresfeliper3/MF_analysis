from load import loader
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.services.OrganismsService import OrganismsService
from src.Biocode.sequences.Sequence import Sequence
from utils.decorators import Timer, DBConnection, TryExcept
from utils.logger import logger


class Analyzer:
    def __init__(self):
        self.organism = ""

    @Timer
    def _load_organism(self, organism_name, gcf, amount_chromosomes):
        logger.info(f"Loading organism {organism_name} - {gcf} - {amount_chromosomes} chromosomes - to database")
        organism_service = OrganismsService()
        organism_service.insert(record=(organism_name, gcf, amount_chromosomes))

    @DBConnection
    @TryExcept
    @Timer
    def analyze_genome_command(self, args):
        save_to_db = False if args.save_to_db == 'false' else True

        if args.name:
            self.organism= args.name
            loader.set_organism(self.organism)
            self._validate_mode_analyzing_genome(args, save_to_db=save_to_db)

        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

    def _validate_mode_analyzing_genome(self, args, save_to_db: bool):
        if args.mode:
            if args.mode == 'whole':
                self._load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                              amount_chromosomes=loader.get_amount_chromosomes())
                self.__whole_MFA_genome(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data(),
                                        save_to_db=save_to_db)
            elif args.mode == 'regions':
                self._load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                              amount_chromosomes=loader.get_amount_chromosomes())
                self.__regions_MFA_genome(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data(),
                                          regions_number=loader.get_regions_number(), save_to_db=save_to_db)
            else:
                raise Exception("Enter a valid mode (whole or regions)")
        else:
            raise Exception("Enter a valid mode (whole or regions)")

    def __whole_MFA_genome(self, organism_name, gcf, data, save_to_db):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db)
        #genome_manager.graph_linear_fit()
        #genome_manager.generate_df_results()

    def __regions_MFA_genome(self, organism_name, gcf, data, regions_number, save_to_db):
        region_genome_manager = RegionGenomeManager(genome_data=data, regions_number=regions_number,
                                                    organism_name=organism_name)
        region_genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db)

    @DBConnection
    @TryExcept
    @Timer
    def analyze_sequence_command(self, args):
        if args.name:
            self.organism= args.name
            loader.set_organism(self.organism)
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

        if args.path:
            sequence = Sequence(sequence=loader.read_fasta_sequence(file_path=args.path),
                                name=loader.extract_file_name(file_path=args.path),
                                organism_name=loader.get_organism_name(),
                                refseq_accession_number=loader.extract_refseq_accession_number(args.path))
            save_to_db = False if args.save_to_db == 'false' else True
            self._validate_mode_analyzing_sequence(args, sequence, save_to_db=save_to_db)
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _validate_mode_analyzing_sequence(self, args, sequence: Sequence, save_to_db: bool):
        if args.mode:
            if args.mode == 'whole':
                self._load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                              amount_chromosomes=loader.get_amount_chromosomes())
                self.__whole_MFA_sequence(gcf=loader.get_gcf(), sequence=sequence, save_to_db=save_to_db)
            elif args.mode == 'regions':
                self._load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                              amount_chromosomes=loader.get_amount_chromosomes())
                self.__regions_MFA_sequence(gcf=loader.get_gcf(), sequence=sequence,
                                     regions_number=loader.get_regions_number(), save_to_db=save_to_db)
            else:
                raise Exception("Enter a valid mode (whole or regions)")
        else:
            raise Exception("Enter a valid mode (whole or regions)")

    def __whole_MFA_sequence(self, gcf, sequence, save_to_db):
        sequence_manager = SequenceManager(sequence=sequence)
        sequence_manager.calculate_multifractal_analysis_values(gcf)

        if save_to_db:
            sequence_manager.save_results_to_db_during_execution(GCF=gcf)

    def __regions_MFA_sequence(self, gcf, sequence, regions_number, save_to_db):
        region_sequence_manager = RegionSequenceManager(sequence=sequence, regions_number=regions_number)
        region_sequence_manager.calculate_multifractal_analysis_values(gcf)

        if save_to_db:
            region_sequence_manager.save_results_to_db_during_execution(GCF=gcf)

    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_genome_command(self, args):
        save_to_db = False if args.save_to_db == 'false' else True

        if args.method:
            if args.method == 'r':
                self.organism= args.name
                loader.set_organism(self.organism)
                self._load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                              amount_chromosomes=loader.get_amount_chromosomes())
                self._find_kmers_recursively_in_genome(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                                                 data=loader.get_data(), save_to_db=save_to_db)
            elif args.method == 'rm':
                logger.warning("Feature not implemented yet")


    def _find_kmers_recursively_in_genome(self, organism_name, gcf, data, save_to_db):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        genome_manager.find_only_kmers_recursively_and_calculate_multifractal_analysis_values(
            GCF=gcf, save_to_db=save_to_db, method_to_find_it="Recursively")

    @DBConnection
    @TryExcept
    @Timer
    def find_kmers_sequence_command(self, args):
        if args.name:
            self.organism= args.name
            loader.set_organism(self.organism)
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

        if args.path:
            if args.method == 'r':
                sequence = Sequence(sequence=loader.read_fasta_sequence(file_path=args.path),
                                    name=loader.extract_file_name(file_path=args.path),
                                    organism_name=loader.get_organism_name(),
                                    refseq_accession_number=loader.extract_refseq_accession_number(args.path))
                save_to_db = False if args.save_to_db == 'false' else True

                self._load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                              amount_chromosomes=loader.get_amount_chromosomes())
                self._find_kmers_recursively_in_sequence(gcf=loader.get_gcf(), sequence=sequence, save_to_db=save_to_db)
            elif args.method == 'rm':
                logger.warning("Feature not implemented yet")
        else:
            raise Exception("Please provide a .fasta file path relative to command.py file")

    def _find_kmers_recursively_in_sequence(self, gcf, sequence, save_to_db):
        sequence_manager = SequenceManager(sequence=sequence)
        sequence_manager.calculate_multifractal_analysis_values(gcf)
        kmers_list = sequence_manager.find_only_kmers_recursively()
        if save_to_db:
            sequence_manager.save_repeats_found_recursively_to_db(
                kmers_list=kmers_list, GCF=gcf, method_to_find_it="Recursively")



