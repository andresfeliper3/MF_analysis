import argparse
from load import loader
from utils.logger import logger

from analyze import load_organism, whole_MFA, regions_MFA
from graph import load_data_whole, graph_whole, load_data_regions, graph_regions
from download import remove_files, execute_download_command, clean_directory, uncompress_all_files


def main():
    parser = argparse.ArgumentParser(description='Command line utility')
    subparsers = parser.add_subparsers(title='Commands', dest='command')

    analyze_parser = subparsers.add_parser('analyze', help='Analyze command')
    analyze_parser.add_argument('-name', help='Name or GCF for analysis')
    analyze_parser.add_argument('-mode', help='Analysis mode: whole / regions')

    graph_parser = subparsers.add_parser('graph', help='Graph command')
    graph_parser.add_argument('-name', help='Name or GCFfor graphing')
    graph_parser.add_argument('-mode', help='Analysis mode: whole / regions')

    download_parser = subparsers.add_parser('download', help='Download command')
    download_parser.add_argument('-name', help='Name or GCF for downloading')

    args = parser.parse_args()

    if args.command == 'analyze':
        analyze_command(args)
    elif args.command == 'graph':
        graph_command(args)
    elif args.command == 'download':
        download_command(args)


organism = ""


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


def analyze_command(args):
    global organism

    if args.name:
        organism = args.name
        loader.set_organism(organism)
        _validate_mode_analyzing(args)

    else:
        logger.error("Please provide either a -name (lowercase name or GCF).")


def _validate_mode_analyzing(args):
    if args.mode:
        if args.mode == 'whole':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            whole_MFA(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data())
        elif args.mode == 'regions':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            regions_MFA(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data(),
                        regions_number=loader.get_regions_number())
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
