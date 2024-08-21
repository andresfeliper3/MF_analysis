import argparse
import traceback
from load import loader
from utils.logger import logger

from analyze import load_organism, whole_MFA_genome, regions_MFA_genome, whole_MFA_sequence, regions_MFA_sequence, \
    find_kmers_recursively_in_genome, find_kmers_recursively_in_sequence
from graph import load_data_whole, graph_whole, load_data_regions, graph_regions, graph_rm_results_from_file, \
    graph_rm_results_from_database, graph_recursive_from_database, graph_recursive_genome_from_database, \
    graph_rm_results_from_files_in_folder, graph_rm_results_of_genome_from_database, graph_gtf_from_file
from download import remove_files, execute_download_command, clean_directory, uncompress_all_files
from repeats import load_RM_repeats_from_file, load_RM_repeats_from_folder
from genes import load_genes_from_file

from src.Biocode.sequences.Sequence import Sequence

def main():
    parser = argparse.ArgumentParser(description='Command line utility')
    subparsers = parser.add_subparsers(title='Commands', dest='command')

    analyze_genome_parser = subparsers.add_parser('analyze_genome', help='Analyze command for a whole genome')
    analyze_genome_parser.add_argument('-name', help='Name or GCF for analysis')
    analyze_genome_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    analyze_genome_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')

    analyze_sequence_parser = subparsers.add_parser('analyze_sequence', help='Analyze command for a sequence')
    analyze_sequence_parser.add_argument('-path', help='Path of the .fasta sequence file relative to command.py file')
    analyze_sequence_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    analyze_sequence_parser.add_argument('-name', help='Name or GCF for analysis')
    analyze_sequence_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')

    kmers_finder_genome_parser = subparsers.add_parser('find_kmers_genome', help='Find kmers (repeats of k nucleotides)')
    kmers_finder_genome_parser.add_argument('-method', help='Define the method to find the kmers (recursively or repeatmasker)')
    #kmers_finder_genome_parseradd_argument('-mode', help='Analysis mode: whole / regions')
    kmers_finder_genome_parser.add_argument('-name', help='Name or GCF for analysis')
    kmers_finder_genome_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')


    kmers_finder_sequence_parser = subparsers.add_parser('find_kmers_sequence', help='Find kmers (repeats of k nucleotides)')
    kmers_finder_sequence_parser.add_argument('-path', help='Path of the .fasta sequence file relative to command.py file')
    kmers_finder_sequence_parser.add_argument('-method', help='Define the method to find the kmers (recursively or repeatmasker)')
    # kmers_finder_sequence_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    kmers_finder_sequence_parser.add_argument('-name', help='Name or GCF for analysis')
    kmers_finder_sequence_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')


    graph_parser = subparsers.add_parser('graph', help='Graph command')
    graph_parser.add_argument('-name', help='Name or GCF for graphing')
    graph_parser.add_argument('-mode', help='Analysis mode: whole / regions')

    graph_rm_file_parser = subparsers.add_parser('graph_rm_file_sequence', help='Graph RepeatMasker results')
    graph_rm_file_parser.add_argument('-path', help="Enter the path of a RepeatMasker results file")
    graph_rm_file_parser.add_argument('-partitions', help="Enter the number of partitions to use to divide the sequence and "
                                                     "merge the repeats")
    graph_rm_file_parser.add_argument('-regions', help="Enter the amount of regions to separate the graph using vertical lines")
    graph_rm_file_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_rm_file_parser.add_argument('--save', choices=['true', 'false'], default='true', help='Save graphs locally in /out directory')
    graph_rm_file_parser.add_argument('-name', help="Enter the scientific name of the organism to use it as a folder name")

    graph_rm_database_parser = subparsers.add_parser('graph_rm_database_sequence', help='Graph RepeatMasker results')
    graph_rm_database_parser.add_argument('-ran', help="Enter the refseq accession number of the sequence/chromosome")
    graph_rm_database_parser.add_argument('-partitions',
                                      help="Enter the number of partitions to use to divide the sequence and "
                                           "merge the repeats")
    graph_rm_database_parser.add_argument('-regions',
                                      help="Enter the amount of regions to separate the graph using vertical lines")
    graph_rm_database_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_rm_database_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                      help='Save graphs locally in /out directory')
    graph_rm_database_parser.add_argument('-name',
                                      help="Enter the scientific name of the organism to use it as a folder name")


    graph_rm_file_genome_parser = subparsers.add_parser('graph_rm_file_genome', help='Graph RepeatMasker results from files in a folder')
    graph_rm_file_genome_parser.add_argument('-path', help="Enter the path of a folder that contains RM results files")
    graph_rm_file_genome_parser.add_argument('-partitions',
                                      help="Enter the number of partitions to use to divide the sequences and "
                                           "merge the repeats")
    graph_rm_file_genome_parser.add_argument('-regions',
                                      help="Enter the amount of regions to separate the graph using vertical lines")
    graph_rm_file_genome_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_rm_file_genome_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                      help='Save graphs locally in /out directory')
    graph_rm_file_genome_parser.add_argument('-name',
                                      help="Enter the scientific name of the organism to use it as a folder name")

    graph_rm_database_genome_parser = subparsers.add_parser('graph_rm_database_genome',
                                                        help='Graph RepeatMasker results from database')
    graph_rm_database_genome_parser.add_argument('-gcf', help="Enter the GCF of the organism genome")
    graph_rm_database_genome_parser.add_argument('-partitions',
                                             help="Enter the number of partitions to use to divide the sequences and "
                                                  "merge the repeats")
    graph_rm_database_genome_parser.add_argument('-regions',
                                             help="Enter the amount of regions to separate the graph using vertical lines")
    graph_rm_database_genome_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_rm_database_genome_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                             help='Save graphs locally in /out directory')
    graph_rm_database_genome_parser.add_argument('-name',
                                             help="Enter the scientific name of the organism to use it as a folder name")

    graph_recursive_parser = subparsers.add_parser('graph_recursive', help='Graph results found by Recursive algorithm')
    graph_recursive_parser.add_argument('-ran', help="Enter the refseq accession number of the sequence/chromosome")
    graph_recursive_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                          help='Save graphs locally in /out directory')
    graph_recursive_parser.add_argument('-name',
                                          help="Enter the scientific name of the organism to use it as a folder name")
    graph_recursive_parser.add_argument('-n_max', help="(Optional) Graph only the top n largest values")

    graph_recursive_genome_parser = subparsers.add_parser('graph_recursive_genome', help='Graph results found by '
                                                                                         'Recursive algorithm in an organism')
    graph_recursive_genome_parser.add_argument('-gcf', help="Enter the GCF id of the organism")
    graph_recursive_genome_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                          help='Save graphs locally in /out directory')
    graph_recursive_genome_parser.add_argument('-name',
                                          help="Enter the scientific name of the organism to use it as a folder name")
    graph_recursive_genome_parser.add_argument('-n_max', help="(Optional) Graph only the top n largest values per sequence")


    graph_gtf_file_parser = subparsers.add_parser('graph_gtf_file', help='Graph genes from .gtf file')
    graph_gtf_file_parser.add_argument('-path', help="Enter the path of a genes .gtf file")
    graph_gtf_file_parser.add_argument('-partitions',
                                      help="Enter the number of partitions to use to divide the sequence and "
                                           "merge the repeats")
    graph_gtf_file_parser.add_argument('-regions',
                                      help="Enter the amount of regions to separate the graph using vertical lines")
    graph_gtf_file_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_gtf_file_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                      help='Save graphs locally in /out directory')
    graph_gtf_file_parser.add_argument('-name',
                                      help="Enter the scientific name of the organism to use it as a folder name")

    download_parser = subparsers.add_parser('download', help='Download command: it downloads the chromosomes files form the link'
                                                             'specified in the sequences.yaml file.')
    download_parser.add_argument('-name', help='Name or GCF for downloading')
    download_parser.add_argument('--gff', choices=['true', 'false'], default='true',
                                      help='Save the genes .gff file from the link specified in the sequence.yaml file')

    load_RM_repeats_parser = subparsers.add_parser('load_RM_repeats', help='Load repeats from results file generated '
                                            'by RepeatMasker')
    load_RM_repeats_parser.add_argument('-path', help='Add the relative path of the .out file. For example: '
                                                      'resources/RM_resources/caenorhabditis_elegans/*_chromosome_I.fasta.out')

    load_RM_repeats_folder_parser = subparsers.add_parser('load_RM_repeats_folder', help='Load repeats from folder of RM results to database')
    load_RM_repeats_folder_parser.add_argument('-path', help='Add the relative path of the folder of RM results. For example: '
                                                      'resources/RM_resources/caenorhabditis_elegans')

    load_genes_parser = subparsers.add_parser('load_genes', help='Load genes file to database using a file path')
    load_genes_parser.add_argument('-path', help='Add the relative path of the genes file. For example: '
                                                      'resources/dna_sequences/caenorhabditis_elegans/gtf/file.gtf')

    args = parser.parse_args()

    if args.command == 'analyze_genome':
        analyze_genome_command(args)
    elif args.command == 'analyze_sequence':
        analyze_sequence_command(args)
    elif args.command == 'find_kmers_genome':
        find_kmers_genome_command(args)
    elif args.command == 'find_kmers_sequence':
        find_kmers_sequence_command(args)
    elif args.command == 'graph':
        graph_command(args)
    elif args.command == 'graph_rm_file_sequence':
        graph_rm_file_command(args)
    elif args.command == 'graph_rm_database_sequence':
        graph_rm_database_command(args)
    elif args.command == 'graph_rm_file_genome':
        graph_rm_file_genome_command(args)
    elif args.command == 'graph_rm_database_genome':
        graph_rm_database_genome_command(args)
    elif args.command == 'graph_recursive':
        graph_recursive_command(args)
    elif args.command == 'graph_recursive_genome':
        graph_recursive_genome_command(args)
    elif args.command == 'graph_gtf_file':
        graph_gtf_file(args)
    elif args.command == 'download':
        download_command(args)
    elif args.command == 'load_RM_repeats':
        load_RM_repeats(args)
    elif args.command == 'load_RM_repeats_folder':
        load_RM_repeats_folder(args)
    elif args.command == 'load_genes':
        load_genes(args)




organism = ""

def find_kmers_genome_command(args):
    global organism

    save_to_db = False if args.save_to_db == 'false' else True

    if args.method:
        if args.method == 'r':
            organism = args.name
            loader.set_organism(organism)
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            find_kmers_recursively_in_genome(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                                             data=loader.get_data(),save_to_db=save_to_db)
        elif args.method == 'rm':
            logger.warning("Feature not implemented yet")


def find_kmers_sequence_command(args):
    global organism

    if args.name:
        organism = args.name
        loader.set_organism(organism)
    else:
        logger.error("Please provide either a -name (lowercase name or GCF).")

    if args.path:
        if args.method == 'r':
            sequence = Sequence(sequence=loader.read_fasta_sequence(file_path=args.path),
                                name=loader.extract_file_name(file_path=args.path),
                                organism_name=loader.get_organism_name(),
                                refseq_accession_number=loader.extract_refseq_accession_number(args.path))
            save_to_db = False if args.save_to_db == 'false' else True

            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            find_kmers_recursively_in_sequence(organism_name=loader.get_organism_name(),
                               sequence_name=loader.extract_file_name(file_path=args.path),
                               gcf=loader.get_gcf(), sequence=sequence, save_to_db=save_to_db)
        elif args.method == 'rm':
            logger.warning("Feature not implemented yet")
    else:
        logger.error("Please provide a .fasta file path relative to command.py file")

def download_command(args):
    global organism

    if args.name:
        organism = args.name
        loader.set_organism(organism)
        remove_files(folder=loader.get_organism_folder())
        execute_download_command(folder=loader.get_organism_folder(), download_url=loader.get_download_url(),
                                 suffix=".gz")
        uncompress_all_files(folder=loader.get_organism_folder())
        if bool(args.gff):
            execute_download_command(folder=loader.get_organism_gtf_subfolder(),
                                     download_url=loader.get_download_gff_url(), suffix=".gtf.gz")
            uncompress_all_files(folder=loader.get_organism_gtf_subfolder())
    else:
        logger.error("Please provide either a -name (lowercase name or GCF).")


def analyze_genome_command(args):
    global organism

    save_to_db = False if args.save_to_db == 'false' else True

    if args.name:
        organism = args.name
        loader.set_organism(organism)
        _validate_mode_analyzing_genome(args, save_to_db=save_to_db)

    else:
        logger.error("Please provide either a -name (lowercase name or GCF).")

def _validate_mode_analyzing_genome(args, save_to_db: bool):
    if args.mode:
        if args.mode == 'whole':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            whole_MFA_genome(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data(),
                             save_to_db=save_to_db)
        elif args.mode == 'regions':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            regions_MFA_genome(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data(),
                               regions_number=loader.get_regions_number(), save_to_db=save_to_db)
        else:
            logger.error("Enter a valid mode (whole or regions)")
    else:
        logger.error("Enter a valid mode (whole or regions)")


def analyze_sequence_command(args):
    global organism

    if args.name:
        organism = args.name
        loader.set_organism(organism)
    else:
        logger.error("Please provide either a -name (lowercase name or GCF).")

    if args.path:
        sequence = Sequence(sequence=loader.read_fasta_sequence(file_path=args.path),
                            name=loader.extract_file_name(file_path=args.path),
                            organism_name=loader.get_organism_name(),
                            refseq_accession_number=loader.extract_refseq_accession_number(args.path))
        save_to_db = False if args.save_to_db == 'false' else True
        _validate_mode_analyzing_sequence(args, sequence, save_to_db=save_to_db)
    else:
        logger.error("Please provide a .fasta file path relative to command.py file")


def _validate_mode_analyzing_sequence(args, sequence: Sequence, save_to_db: bool):
    if args.mode:
        if args.mode == 'whole':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            whole_MFA_sequence(gcf=loader.get_gcf(), sequence=sequence, save_to_db=save_to_db)
        elif args.mode == 'regions':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            regions_MFA_sequence(gcf=loader.get_gcf(), sequence=sequence,
                                 regions_number=loader.get_regions_number(), save_to_db=save_to_db)
        else:
            logger.error("Enter a valid mode (whole or regions)")
    else:
        logger.error("Enter a valid mode (whole or regions)")

def graph_command(args):
    global organism

    if args.name:
        organism = args.name
        loader.set_organism(organism)
        _validate_mode_graphing(args)
    else:
        logger.error("Please provide either -id or -name.")


def _validate_mode_graphing(args):
    if args.mode:
        if args.mode == 'whole':
            df = load_data_whole(gcf=loader.get_gcf())
            graph_whole(dataframe=df, organism_name=loader.get_organism_name(), data=loader.get_data())
        elif args.mode == 'regions':
            df = load_data_regions(gcf=loader.get_gcf())
            graph_regions(dataframe=df, organism_name=loader.get_organism_name(), data=loader.get_data(),
                          regions_number=loader.get_regions_number())
        else:
            logger.error("Enter a valid mode (whole or regions)")
    else:
        logger.error("Enter a valid mode (whole or regions)")


def load_RM_repeats(args):
    try:
        load_RM_repeats_from_file(args.path)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()


def load_RM_repeats_folder(args):
    try:
        load_RM_repeats_from_folder(args.path)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()


def graph_rm_file_command(args):
    try:
        graph_rm_results_from_file(path=args.path, partitions=args.partitions, regions=args.regions,
                                   plot_type=args.plot_type, save=args.save, name=args.name)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

def graph_rm_database_command(args):
    try:
        graph_rm_results_from_database(refseq_accession_number=args.ran, partitions=args.partitions,
                                   regions=args.regions, plot_type=args.plot_type, save=args.save, name=args.name)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()


def graph_rm_file_genome_command(args):
    try:
        graph_rm_results_from_files_in_folder(directory_path=args.path, partitions=args.partitions, regions=args.regions,
                                              plot_type=args.plot_type, save=args.save, name=args.name)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

def graph_rm_database_genome_command(args):
    try:
        graph_rm_results_of_genome_from_database(GCF=args.gcf, partitions=args.partitions, regions=args.regions,
                                                 plot_type=args.plot_type, save=args.save, name=args.name)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

def graph_recursive_command(args):
    try:
        graph_recursive_from_database(refseq_accession_number=args.ran, save=args.save, name=args.name, n_max=args.n_max)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

def graph_recursive_genome_command(args):
    try:
        graph_recursive_genome_from_database(GCF=args.gcf, save=args.save, name=args.name, n_max=args.n_max)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

def graph_gtf_file(args):
    try:
        graph_gtf_from_file(path=args.path, partitions=args.partitions, regions=args.regions, plot_type=args.plot_type,
                            save=args.save, name=args.name)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

def  load_genes(args):
    try:
        load_genes_from_file(path=args.path)
    except Exception as e:
        logger.error(e)
        traceback.print_exc()

if __name__ == "__main__":
    main()

