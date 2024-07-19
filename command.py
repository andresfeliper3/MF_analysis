import argparse
from load import loader
from utils.logger import logger

from analyze import load_organism, whole_MFA_genome, regions_MFA_genome, whole_MFA_sequence, regions_MFA_sequence, \
    find_kmers_recursively_in_genome, find_kmers_recursively_in_sequence
from graph import load_data_whole, graph_whole, load_data_regions, graph_regions
from download import remove_files, execute_download_command, clean_directory, uncompress_all_files

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

    download_parser = subparsers.add_parser('download', help='Download command')
    download_parser.add_argument('-name', help='Name or GCF for downloading')

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
    elif args.command == 'download':
        download_command(args)



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
            sequence = Sequence(sequence=loader.read_fasta_sequence(file_path=args.path))
            save_to_db = False if args.save_to_db == 'false' else True

            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            find_kmers_recursively_in_sequence(organism_name=loader.get_organism_name(),
                               sequence_name=loader.extract_sequence_name(file_path=args.path),
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
        logger.debug(loader.get_organism_folder())
        remove_files(organism_folder=loader.get_organism_folder())
        execute_download_command(organism_folder=loader.get_organism_folder(), download_url=loader.get_download_url())
        clean_directory(organism_folder=loader.get_organism_folder())
        uncompress_all_files(organism_folder=loader.get_organism_folder())
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
        sequence = Sequence(sequence=loader.read_fasta_sequence(file_path=args.path))
        save_to_db = False if args.save_to_db == 'false' else True
        _validate_mode_analyzing_sequence(args, sequence, save_to_db=save_to_db)
    else:
        logger.error("Please provide a .fasta file path relative to command.py file")


def _validate_mode_analyzing_sequence(args, sequence: Sequence, save_to_db: bool):
    if args.mode:
        if args.mode == 'whole':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            whole_MFA_sequence(organism_name=loader.get_organism_name(),
                               sequence_name=loader.extract_sequence_name(file_path=args.path),
                               gcf=loader.get_gcf(), sequence=sequence, save_to_db=save_to_db)
        elif args.mode == 'regions':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            regions_MFA_sequence(organism_name=loader.get_organism_name(),
                                 sequence_name=loader.extract_sequence_name(file_path=args.path),
                                 gcf=loader.get_gcf(), sequence=sequence,
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


if __name__ == "__main__":
    main()
