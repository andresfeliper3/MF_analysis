import argparse
from load import loader
from utils.logger import logger


def main():
    parser = argparse.ArgumentParser(description='Command line utility')
    subparsers = parser.add_subparsers(title='Commands', dest='command')

    analyze_parser = subparsers.add_parser('analyze', help='Analyze command')
    analyze_parser.add_argument('-id', help='ID for analysis')
    analyze_parser.add_argument('-name', help='Name for analysis')

    graph_parser = subparsers.add_parser('graph', help='Graph command')
    graph_parser.add_argument('-id', help='ID for graphing')
    graph_parser.add_argument('-name', help='Name for graphing')

    args = parser.parse_args()

    if args.command == 'analyze':
        analyze_command(args)
    elif args.command == 'graph':
        graph_command(args)


organism = ""


def analyze_command(args):
    global organism
    if args.id:
        organism = args.id
        loader.set_organism(organism)
        logger.debug(loader.get_data())

    elif args.name:
        organism = args.name
        loader.set_organism(organism)
        logger.debug(loader.get_data())

    else:
        print("Please provide either -id or -name.")


def graph_command(args):
    if args.id:
        print("Graphing by ID:", args.id)
    elif args.name:
        print("Graphing by name:", args.name)
    else:
        print("Please provide either -id or -name.")


if __name__ == "__main__":
    main()
