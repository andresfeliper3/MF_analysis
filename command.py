import argparse
from load import loader
from utils.logger import logger

from analyze import load_organism, whole_MFA, regions_MFA


def main():
    parser = argparse.ArgumentParser(description='Command line utility')
    subparsers = parser.add_subparsers(title='Commands', dest='command')

    analyze_parser = subparsers.add_parser('analyze', help='Analyze command')
    analyze_parser.add_argument('-name', help='Name for analysis')
    analyze_parser.add_argument('-mode', help='Analysis mode: whole / regions')

    graph_parser = subparsers.add_parser('graph', help='Graph command')
    graph_parser.add_argument('-name', help='Name for graphing')

    args = parser.parse_args()

    if args.command == 'analyze':
        analyze_command(args)
    elif args.command == 'graph':
        graph_command(args)


organism = ""


def analyze_command(args):
    global organism

    if args.name:
        organism = args.name
        loader.set_organism(organism)
        _validate_mode(args)

    else:
        logger.error("Please provide either a -name (lowercase name or GCF).")


def _validate_mode(args):
    if args.mode:
        if args.mode == 'whole':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            whole_MFA(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data())
        elif arg.mode == 'regions':
            load_organism(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(),
                          amount_chromosomes=loader.get_amount_chromosomes())
            regions_MFA(organism_name=loader.get_organism_name(), gcf=loader.get_gcf(), data=loader.get_data(),
                        regions_number=loader.get_regions_number())
        else:
            logger.error("Enter a valid mode (whole or regions)")
    else:
        logger.error("Enter a valid mode (whole or regions)")


def graph_command(args):
    if args.id:
        print("Graphing by ID:", args.id)
    elif args.name:
        print("Graphing by name:", args.name)
    else:
        print("Please provide either -id or -name.")


if __name__ == "__main__":
    main()
