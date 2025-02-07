from Loader import Loader
from src.Biocode.managers.GenomeManager import GenomeManager
from src.Biocode.managers.RegionGenomeManager import RegionGenomeManager
from src.Biocode.managers.RegionSequenceManager import RegionSequenceManager
from src.Biocode.managers.SequenceManager import SequenceManager
from src.Biocode.sequences.Sequence import Sequence
from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.services.OrganismsService import OrganismsService
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.logger import logger


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
            self.organism = args.name
            self.loader.set_organism(self.organism)
            dir_name = args.name.replace(" ", "_")
            self._validate_mode_analyzing_genome(args, save_to_db=save_to_db, name=dir_name,
                                                 only_cgr=bool(args.only_cgr))
        else:
            raise Exception("Please provide either a -name (lowercase name or GCF).")

    def _validate_mode_analyzing_genome(self, args, save_to_db: bool, name: str, only_cgr: bool):
        if args.mode:
            if args.mode == 'whole':
                self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                              amount_chromosomes=self.loader.get_amount_chromosomes())
                self.__whole_MFA_genome(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(), data=self.loader.get_data(),
                                        save_to_db=save_to_db, name=name, only_cgr=only_cgr)
            elif args.mode == 'regions':
                self._load_organism(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(),
                              amount_chromosomes=self.loader.get_amount_chromosomes())
                self.__regions_MFA_genome(organism_name=self.loader.get_organism_name(), gcf=self.loader.get_gcf(), data=self.loader.get_data(),
                                          regions_number=args.regions_number, window_length=args.window_length,
                                          save_to_db=save_to_db, name=name, only_cgr=only_cgr)
            else:
                raise Exception("Enter a valid mode (whole or regions)")
        else:
            raise Exception("Enter a valid mode (whole or regions)")

    def __whole_MFA_genome(self, organism_name, gcf, data, save_to_db, name, only_cgr):
        genome_manager = GenomeManager(genome_data=data, organism_name=organism_name)
        genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db, name=name,
                                                              only_cgr=only_cgr)
        #genome_manager.graph_linear_fit()
        #genome_manager.generate_df_results()

    def __regions_MFA_genome(self, organism_name, gcf, data, regions_number, window_length, save_to_db, name, only_cgr):
        window_length = int(window_length) if window_length else None
        regions_number = int(regions_number) if regions_number else None
        region_genome_manager = RegionGenomeManager(genome_data=data, regions_number=regions_number,
                                                    window_length=window_length, organism_name=organism_name)
        region_genome_manager.calculate_multifractal_analysis_values(GCF=gcf, save_to_db=save_to_db, name=name,
                                                                     only_cgr=only_cgr)

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



