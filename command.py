import argparse

from Analyzer import Analyzer
from Downloader import Downloader
from GenesLoader import GenesLoader
from Grapher import Grapher
from RepeatsLoader import RepeatsLoader


def main():
    parser = argparse.ArgumentParser(description='Command line utility')
    subparsers = parser.add_subparsers(title='Commands', dest='command')

    analyze_genome_parser = subparsers.add_parser('analyze_genome', help='Analyze command for a whole genome')
    analyze_genome_parser.add_argument('-name', help='Name or GCF for analysis')
    analyze_genome_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    analyze_genome_parser.add_argument('-regions_number', help='Enter the amount of regions')
    analyze_genome_parser.add_argument('-window_length', help='Enter the length of the windows')
    analyze_genome_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')

    analyze_sequence_parser = subparsers.add_parser('analyze_sequence', help='Analyze command for a sequence')
    analyze_sequence_parser.add_argument('-path', help='Path of the .fasta sequence file relative to command.py file')
    analyze_sequence_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    analyze_sequence_parser.add_argument('-regions_number', help='Enter the amount of regions')
    analyze_sequence_parser.add_argument('-window_length', help='Enter the length of the windows')
    analyze_sequence_parser.add_argument('-name', help='Name or GCF for analysis')
    analyze_sequence_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')

    kmers_finder_genome_parser = subparsers.add_parser('find_kmers_genome', help='Find kmers (repeats of k nucleotides)')
    kmers_finder_genome_parser.add_argument('-method', help='Define the method to find the kmers (recursively or repeatmasker)')
    #kmers_finder_genome_parseradd_argument('-mode', help='Analysis mode: whole / regions')
    kmers_finder_genome_parser.add_argument('-k_range', help="Add a k_range to find kmers of k values between the interval. For example: (4, 8).")
    kmers_finder_genome_parser.add_argument('-name', help='Name or GCF for analysis')
    kmers_finder_genome_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')
    kmers_finder_genome_parser.add_argument('-dir', help='Directory to save results')
    kmers_finder_genome_parser.add_argument('-window_length', help='Window/partition length to divide the chromosome')
    kmers_finder_genome_parser.add_argument('-graph_from_file', choices=['true', 'false'], default='true',
                                              help='Generate graphs directly from file without database')

    kmers_finder_sequence_parser = subparsers.add_parser('find_kmers_sequence', help='Find kmers (repeats of k nucleotides)')
    kmers_finder_sequence_parser.add_argument('-path', help='Path of the .fasta sequence file relative to command.py file')
    kmers_finder_sequence_parser.add_argument('-method', help='Define the method to find the kmers (recursively or linearly) (-r or -l)')
    # kmers_finder_sequence_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    kmers_finder_sequence_parser.add_argument('-k_range', help="Add a k_range to find kmers of k values between the interval. For example: (4, 8).")
    kmers_finder_sequence_parser.add_argument('-name', help='Name or GCF for analysis')
    kmers_finder_sequence_parser.add_argument('-dir', help='Directory to save results')
    kmers_finder_sequence_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true', help='Save results to the database')
    kmers_finder_sequence_parser.add_argument('-window_length', help='Window/partition length to divide the chromosome')
    kmers_finder_sequence_parser.add_argument('-graph_from_file', choices=['true', 'false'], default='true',
                                                      help='Generate graphs directly from file without database')

    kmers_finder_linearly_genes_sequence_parser = subparsers.add_parser('find_kmers_linearly_genes_sequence', help='Find kmers using linear method(repeats of k nucleotides) only in the genes')
    kmers_finder_linearly_genes_sequence_parser.add_argument('-path',
                                              help='Path of the .fasta sequence file relative to command.py file')
    kmers_finder_linearly_genes_sequence_parser.add_argument('-name', help='Name or GCF for analysis')
    kmers_finder_linearly_genes_sequence_parser.add_argument('-dir', help='Directory to save results')
    kmers_finder_linearly_genes_sequence_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true',
                                              help='Save results to the database')
    kmers_finder_linearly_genes_sequence_parser.add_argument('-window_length', help='Window/partition length to divide the chromosome')
    kmers_finder_linearly_genes_sequence_parser.add_argument('-graph_from_file', choices=['true', 'false'], default='true',
                                              help='Graph frequency of kmers across sequence only in genes')
    kmers_finder_linearly_genes_sequence_parser.add_argument('-size', help='Size of kmers/repeats that are going to be searched in genes')

    kmers_finder_linearly_genes_genome_parser = subparsers.add_parser('find_kmers_linearly_genes_genome',
                                                                        help='Find kmers using linear method(repeats of k nucleotides) only in the genes')
    kmers_finder_linearly_genes_genome_parser.add_argument('-name', help='Name or GCF for analysis')
    kmers_finder_linearly_genes_genome_parser.add_argument('-dir', help='Directory to save results')
    kmers_finder_linearly_genes_genome_parser.add_argument('--save-to-db', choices=['true', 'false'], default='true',
                                                             help='Save results to the database')
    kmers_finder_linearly_genes_genome_parser.add_argument('-window_length',
                                                             help='Window/partition length to divide the chromosome')
    kmers_finder_linearly_genes_genome_parser.add_argument('-graph_from_file', choices=['true', 'false'],
                                                             default='true',
                                                             help='Graph frequency of kmers across sequence only in genes')
    kmers_finder_linearly_genes_genome_parser.add_argument('-size', help='Size of kmers/repeats that are going to be searched in genes')


    graph_parser = subparsers.add_parser('graph', help='Graph command')
    graph_parser.add_argument('-name', help='Name or GCF for graphing')
    graph_parser.add_argument('-mode', help='Analysis mode: whole / regions')
    graph_parser.add_argument('-regions_number', help='Amount of regions in which the sequence will be divided.')
    graph_parser.add_argument('-window_length', help='Size (amount of bps) in which the windows/regions will be divided.')

    graph_rm_file_parser = subparsers.add_parser('graph_rm_file_sequence', help='Graph RepeatMasker results')
    graph_rm_file_parser.add_argument('-path', help="Enter the path of a RepeatMasker results file")
    graph_rm_file_parser.add_argument('-partitions', help="Enter the number of partitions to use to divide the sequence and "
                                                     "merge the repeats")
    graph_rm_file_parser.add_argument('-regions', help="Enter the amount of regions to separate the graph using vertical lines")
    graph_rm_file_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_rm_file_parser.add_argument('--save', choices=['true', 'false'], default='true', help='Save graphs locally in /out directory')
    graph_rm_file_parser.add_argument('-dir', help="Enter the scientific name of the organism to use it as a folder name")

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
    graph_rm_database_parser.add_argument('-dir',
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
    graph_rm_file_genome_parser.add_argument('-dir',
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
    graph_rm_database_genome_parser.add_argument('-dir',
                                             help="Enter the scientific name of the organism to use it as a folder name")

    graph_genome_repeats_from_file = subparsers.add_parser('graph_genome_repeats_from_file', help='Graph genome repeats results from'
                                                                                        'a single file (such as Plantrep.cn results)')
    graph_genome_repeats_from_file.add_argument('-path', help="Enter the path of a genome results file")
    graph_genome_repeats_from_file.add_argument('-partitions',
                                      help="Enter the number of partitions to use to divide the sequence and "
                                           "merge the repeats")
    graph_genome_repeats_from_file.add_argument('-regions',
                                      help="Enter the amount of regions to separate the graph using vertical lines")
    graph_genome_repeats_from_file.add_argument('-plot_type', help="Plot type: line or bar")
    graph_genome_repeats_from_file.add_argument('--save', choices=['true', 'false'], default='true',
                                      help='Save graphs locally in /out directory')
    graph_genome_repeats_from_file.add_argument('-dir',
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
    graph_recursive_genome_parser.add_argument('-dir',
                                          help="Enter the scientific name of the organism to use it as a folder name")
    graph_recursive_genome_parser.add_argument('-n_max', help="(Optional) Graph only the top n largest values per sequence")

    graph_linear_sequence_repeats_parser = subparsers.add_parser('graph_linear_repeats_sequence', help='Graph results found by Linear algorithm')
    graph_linear_sequence_repeats_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                        help='Save graphs locally in /out directory')
    graph_linear_sequence_repeats_parser.add_argument('-name',
                                        help="Enter the scientific name of the organism")
    graph_linear_sequence_repeats_parser.add_argument('-path', help="Enter the path of the sequence")
    graph_linear_sequence_repeats_parser.add_argument('-dir', help="Enter the directory name (organism name) were the results will be saved")
    graph_linear_sequence_repeats_parser.add_argument('-k_range', help="Add a k_range to graph kmers of k values between the interval. For example: (4, 8).")


    graph_linear_genome_repeats_parser = subparsers.add_parser('graph_linear_repeats_genome',
                                                        help='Graph results found by Linear algorithm')
    graph_linear_genome_repeats_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                             help='Save graphs locally in /out directory')
    graph_linear_genome_repeats_parser.add_argument('-name',
                                             help="Enter the scientific name of the organism")
    graph_linear_genome_repeats_parser.add_argument('-gcf', help="Enter the GCF of the organism genome")
    graph_linear_genome_repeats_parser.add_argument('-dir',
                                             help="Enter the directory name (organism name) were the results will be saved")
    graph_linear_genome_repeats_parser.add_argument('-k_range', help="Add a k_range to graph kmers of k values between the interval. For example: (4, 8).")



    graph_linear_in_genes_sequence_repeats_parser = subparsers.add_parser('graph_linear_in_genes_repeats_sequence',
                                                                 help='Graph results found by Linear algorithm only in genes')
    graph_linear_in_genes_sequence_repeats_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                                      help='Save graphs locally in /out directory')
    graph_linear_in_genes_sequence_repeats_parser.add_argument('-name',
                                                      help="Enter the scientific name of the organism")
    graph_linear_in_genes_sequence_repeats_parser.add_argument('-path', help="Enter the path of the sequence")
    graph_linear_in_genes_sequence_repeats_parser.add_argument('-dir',
                                                      help="Enter the directory name (organism name) were the results will be saved")
    graph_linear_in_genes_sequence_repeats_parser.add_argument('-k_range', help="Add a k_range to graph kmers of k values between the interval. For example: (4, 8).")


    graph_linear_in_genes_genome_repeats_parser = subparsers.add_parser('graph_linear_in_genes_repeats_genome',
                                                               help='Graph results found by Linear algorithm')
    graph_linear_in_genes_genome_repeats_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                                    help='Save graphs locally in /out directory')
    graph_linear_in_genes_genome_repeats_parser.add_argument('-name',
                                                    help="Enter the scientific name of the organism")
    graph_linear_in_genes_genome_repeats_parser.add_argument('-gcf', help="Enter the GCF of the organism genome")
    graph_linear_in_genes_genome_repeats_parser.add_argument('-dir',
                                                    help="Enter the directory name (organism name) were the results will be saved")
    graph_linear_in_genes_genome_repeats_parser.add_argument('-k_range',
                                                    help="Add a k_range to graph kmers of k values between the interval. For example: (4, 8).")


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
    graph_gtf_file_parser.add_argument('-dir',
                                      help="Enter the scientific name of the organism to use it as a folder name")

    graph_gtf_file_parser = subparsers.add_parser('graph_gtf_database', help='Graph genes from .gtf file')
    graph_gtf_file_parser.add_argument('-gcf', help="Enter the GCF of the organism/genome")
    graph_gtf_file_parser.add_argument('-ran', help="Enter the refseq accession number of the sequence/chromosome")
    graph_gtf_file_parser.add_argument('-partitions',
                                       help="Enter the number of partitions to use to divide the sequence and "
                                            "merge the repeats")
    graph_gtf_file_parser.add_argument('-regions',
                                       help="Enter the amount of regions to separate the graph using vertical lines")
    graph_gtf_file_parser.add_argument('-plot_type', help="Plot type: line or bar")
    graph_gtf_file_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                       help='Save graphs locally in /out directory')
    graph_gtf_file_parser.add_argument('-dir',
                                       help="Enter the scientific name of the organism to use it as a folder name")
    
    graph_linear_regression_sequence_parser = subparsers.add_parser('graph_linear_regression_sequence', help='Graph linear regression of DDq vs Repeats')
    graph_linear_regression_sequence_parser.add_argument('-path', help='Path of the chromosome fasta file')
    graph_linear_regression_sequence_parser.add_argument('-k_range', help='Range of values for k in kmers. Ex: "(4,12)".')
    graph_linear_regression_sequence_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                       help='Save graphs locally in /out directory')
    graph_linear_regression_sequence_parser.add_argument('-dir',
                                       help="Enter the scientific name of the organism to use it as a folder name")
    graph_linear_regression_sequence_parser.add_argument('-name', help="Enter the scientific name of the organism")


    graph_linear_regression_genome_parser = subparsers.add_parser('graph_linear_regression_genome',
                                                    help='Graph linear regression of DDq vs Repeats')
    graph_linear_regression_genome_parser.add_argument('-gcf', help='GCF of the organism')
    graph_linear_regression_genome_parser.add_argument('-k_range', help='Range of values for k in kmers. Ex: "(4,12)".')
    graph_linear_regression_genome_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                         help='Save graphs locally in /out directory')
    graph_linear_regression_genome_parser.add_argument('-dir',
                                         help="Enter the scientific name of the organism to use it as a folder name")
    graph_linear_regression_genome_parser.add_argument('-name', help="Enter the scientific name of the organism")


    graph_categories_repeats_heatmap_sequence_parser = subparsers.add_parser('graph_categories_repeats_heatmap_sequence',
                                                                  help='Graph KEGG categories vs repeats heatmap')

    graph_categories_repeats_heatmap_sequence_parser.add_argument('-path', help='Path of the chromosome fasta file')
    graph_categories_repeats_heatmap_sequence_parser.add_argument('-size',
                                                         help='Size value of the kmers/repeats')
    graph_categories_repeats_heatmap_sequence_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                                         help='Save graphs locally in /out directory')
    graph_categories_repeats_heatmap_sequence_parser.add_argument('-dir',
                                                         help="Enter the scientific name of the organism to use it as a folder name")
    graph_categories_repeats_heatmap_sequence_parser.add_argument('-name', help="Enter the scientific name of the organism")
    graph_categories_repeats_heatmap_sequence_parser.add_argument('--tags', choices=['true', 'false'],
                                                                  default='true', help="Include tags/annotations in the heatmap plot")

    graph_categories_repeats_heatmap_genome_parser = subparsers.add_parser('graph_categories_repeats_heatmap_genome',
        help='Graph KEGG categories vs repeats heatmap')

    graph_categories_repeats_heatmap_genome_parser.add_argument('-size',
                                                                  help='Size value of the kmers/repeats')
    graph_categories_repeats_heatmap_genome_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                                                  help='Save graphs locally in /out directory')
    graph_categories_repeats_heatmap_genome_parser.add_argument('-dir',
                                                                  help="Enter the scientific name of the organism to use it as a folder name")
    graph_categories_repeats_heatmap_genome_parser.add_argument('-name',
                                                                  help="Enter the scientific name of the organism")
    graph_categories_repeats_heatmap_genome_parser.add_argument('--tags', choices=['true', 'false'],
                                                                  default='true',
                                                                  help="Include tags/annotations in the heatmap plot")

    graph_compare_parser = subparsers.add_parser('graph_compare',  help='Graph comparing degrees of multifractality of different organisms')
    graph_compare_parser.add_argument('-organisms', nargs='+', help='List of organisms to compare')
    graph_compare_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                                                help='Save graphs locally in /out directory')
    graph_compare_parser.add_argument('-dir',
                                                                help="Enter the name of the genre or family")

    graph_compare_genes_mfa_kmers_parser = subparsers.add_parser('graph_compare_genes_mfa_kmers',
                                                 help='Graph comparing degrees of multifractality of different organisms')
    graph_compare_genes_mfa_kmers_parser.add_argument('-gcf', help='GCF of the organism')
    graph_compare_genes_mfa_kmers_parser.add_argument('-k_range', help='Range of k values for kmers')
    graph_compare_genes_mfa_kmers_parser.add_argument('--save', choices=['true', 'false'], default='true',
                                      help='Save graphs locally in /out directory')
    graph_compare_genes_mfa_kmers_parser.add_argument('-dir',
                                      help="Enter the name of the genre or family")

    download_parser = subparsers.add_parser('download', help='Download command: it downloads the chromosomes files form the link'
                                                             'specified in the sequences.yaml file.')
    download_parser.add_argument('-name', help='Name or GCF for downloading')
    download_parser.add_argument('--gtf', choices=['true', 'false'], default='true',
                                      help='Save the genes .gtf file from the link specified in the sequence.yaml file')

    load_RM_repeats_parser = subparsers.add_parser('load_RM_repeats', help='Load repeats from results file generated '
                                            'by RepeatMasker')
    load_RM_repeats_parser.add_argument('-path', help='Add the relative path of the .out file. For example: '
                                                      'resources/RM_resources/caenorhabditis_elegans/*_chromosome_I.fasta.out')

    load_RM_repeats_folder_parser = subparsers.add_parser('load_RM_repeats_folder', help='Load repeats from folder of RM results to database')
    load_RM_repeats_folder_parser.add_argument('-path', help='Add the relative path of the folder of RM results. For example: '
                                                      'resources/RM_resources/caenorhabditis_elegans')

    load_genome_repeats_file_parser = subparsers.add_parser('load_genome_repeats_file', help='Load repeats of an entire genome'
                                                                                             'from a single results file')
    load_genome_repeats_file_parser.add_argument('-path', help='Add the relative path of the file of genome repeats. For example: '
                                                      'resources/RM_resources/musa_acuminata/157_Musa_acuminata_rm.out')

    load_genes_parser = subparsers.add_parser('load_genes', help='Load genes file to database using a file path')
    load_genes_parser.add_argument('-path', help='Add the relative path of the genes file. For example: '
                                                      'resources/dna_sequences/caenorhabditis_elegans/gtf/file.gtf')

    load_categories_parser = subparsers.add_parser('load_categories',
                                                                        help='Scrape on KEGG website and add categories to corresponding genes')
    load_categories_parser.add_argument('-path',
                                                             help='Path of the .fasta sequence file relative to command.py file')
    load_categories_parser.add_argument('-name', help='Name or GCF for analysis')
    load_categories_parser.add_argument('-size',
                                                    help='Size of kmers/repeats that are going to be searched in genes')
    load_categories_parser.add_argument('-genes_amount', help='The top N amount of genes that will be '
                                                              'considered. Top is based on genes which include more repeats.')

    load_categories_genome_parser = subparsers.add_parser('load_categories_genome',
                                                   help='Scrape on KEGG website and add categories to corresponding genes')
    load_categories_genome_parser.add_argument('-name', help='Name or GCF for analysis')
    load_categories_genome_parser.add_argument('-size',
                                        help='Size of kmers/repeats that are going to be searched in genes')
    load_categories_genome_parser.add_argument('-genes_amount', help='The top N amount of genes that will be '
                                                              'considered. Top is based on genes which include more repeats.')

    args = parser.parse_args()
    analyzer = Analyzer()
    downloader = Downloader()
    genes_loader = GenesLoader()
    repeats_loader = RepeatsLoader()
    grapher = Grapher()

    if args.command == 'analyze_genome':
        analyzer.analyze_genome_command(args)
    elif args.command == 'analyze_sequence':
        analyzer.analyze_sequence_command(args)
    elif args.command == 'find_kmers_genome':
        repeats_loader.find_kmers_genome_command(args)
    elif args.command == 'find_kmers_sequence':
        repeats_loader.find_kmers_sequence_command(args)
    elif args.command == 'find_kmers_linearly_genes_sequence':
        repeats_loader.find_kmers_linearly_genes_sequence_command(filepath=args.path, args=args)
    elif args.command == 'find_kmers_linearly_genes_genome':
        repeats_loader.find_kmers_linearly_genes_genome_command(args)
    elif args.command == 'graph':
        grapher.graph_command(args)
    elif args.command == 'graph_rm_file_sequence':
        grapher.graph_rm_results_from_file(path=args.path, partitions=args.partitions, regions=args.regions,
                                   plot_type=args.plot_type, save=args.save, dir=args.dir)
    elif args.command == 'graph_rm_database_sequence':
        grapher.graph_rm_results_from_database(refseq_accession_number=args.ran, partitions=args.partitions,
                                   regions=args.regions, plot_type=args.plot_type, save=args.save, dir=args.dir)
    elif args.command == 'graph_rm_file_genome':
        grapher.graph_rm_results_from_files_in_folder(directory_path=args.path, partitions=args.partitions, regions=args.regions,
                                             plot_type=args.plot_type, save=args.save, dir=args.dir)
    elif args.command == 'graph_rm_database_genome':
        grapher.graph_rm_results_of_genome_from_database(GCF=args.gcf, partitions=args.partitions, regions=args.regions,
                                                 plot_type=args.plot_type, save=args.save, dir=args.dir)
    elif args.command == 'graph_genome_repeats_from_file':
        grapher.graph_genome_repeats_from_file(path=args.path, dir=args.dir, partitions=args.partitions, regions=args.regions,
                                   plot_type=args.plot_type, save=args.save)
    elif args.command == 'graph_recursive':
        grapher.graph_recursive_from_database(refseq_accession_number=args.ran, save=args.save, name=args.name, n_max=args.n_max)
    elif args.command == 'graph_recursive_genome':
        grapher.graph_recursive_genome_from_database(GCF=args.gcf, save=args.save, dir=args.dir, n_max=args.n_max)
    elif args.command == 'graph_linear_repeats_sequence':
        grapher.graph_linear_repeats_sequence_command(path=args.path, save=args.save, name=args.name, dir=args.dir,
                                                      k_range=args.k_range)
    elif args.command == 'graph_linear_repeats_genome':
        grapher.graph_linear_repeats_genome_command(save=args.save, name=args.name, dir=args.dir, GCF=args.gcf,
                                                    k_range=args.k_range)
    elif args.command == 'graph_linear_in_genes_repeats_sequence':
        grapher.graph_linear_in_genes_repeats_sequence_command(path=args.path, save=args.save, name=args.name,
                                                               dir=args.dir, k_range=args.k_range)
    elif args.command == 'graph_linear_in_genes_repeats_genome':
        grapher.graph_linear_in_genes_repeats_genome_command(GCF=args.gcf, save=args.save, dir=args.dir, name=args.name,
                                                             k_range=args.k_range)
    elif args.command == 'graph_gtf_file':
        grapher.graph_gtf_from_file(path=args.path, partitions=args.partitions, regions=args.regions, plot_type=args.plot_type,
                            save=args.save, dir=args.dir)
    elif args.command == 'graph_gtf_database':
        grapher.graph_gtf_from_database(GCF=args.gcf, refseq_accession_number=args.ran, partitions=args.partitions,
                                regions=args.regions, plot_type=args.plot_type, save=args.save, dir=args.dir)
    elif args.command == 'graph_linear_regression_sequence':
        grapher.graph_linear_regression_sequence_command(k_range=args.k_range, name=args.name, save=args.save,
                                                         dir=args.dir, path=args.path)
    elif args.command == 'graph_linear_regression_genome':
        grapher.graph_linear_regression_genome_command(GCF=args.gcf, k_range=args.k_range, name=args.name, save=args.save,
                                                       dir=args.dir)
    elif args.command == 'graph_categories_repeats_heatmap_sequence':
        grapher.graph_categories_repeats_heatmap_sequence_command(path=args.path, size=args.size, save=args.save, dir=args.dir,
                                                                  name=args.name, tags=args.tags)
    elif args.command == 'graph_categories_repeats_heatmap_genome':
        grapher.graph_categories_repeats_heatmap_genome_command(size=args.size, save=args.save, dir=args.dir,
                                                                  name=args.name, tags=args.tags)
    elif args.command == 'graph_compare':
        grapher.graph_compare_command(organisms=args.organisms, save=args.save, dir=args.dir)
    elif args.command == "graph_compare_genes_mfa_kmers":
        grapher.compare_genes_mfa_kmers_command(GCF=args.gcf, k_range=args.k_range, dir=args.dir, save=args.save)
    elif args.command == 'download':
        downloader.download_command(args)
    elif args.command == 'load_RM_repeats':
        repeats_loader.load_RM_repeats_from_file(args.path)
    elif args.command == 'load_RM_repeats_folder':
        repeats_loader.load_RM_repeats_from_folder(args.path)
    elif args.command == 'load_genome_repeats_file':
       repeats_loader.load_genome_repeats_file(args.path)
    elif args.command == 'load_genes':
        genes_loader.load_genes_from_file(args.path)
    elif args.command == 'load_categories':
        genes_loader.load_categories_from_kegg(filepath=args.path, args=args)
    elif args.command == 'load_categories_genome':
        genes_loader.load_categories_from_kegg_genome(args)

if __name__ == "__main__":
    main()

