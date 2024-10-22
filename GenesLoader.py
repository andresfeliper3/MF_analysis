from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from src.Biocode.services.KeggCategoriesService import KeggCategoriesService
from src.Biocode.services.KeggSubcategoriesService import KeggSubcategoriesService
from src.Biocode.services.GtfGenesKeggCategoriesService import GtfGenesKeggCategoriesService
from src.Biocode.services.GtfGenesKeggSubcategoriesService import GtfGenesKeggSubcategoriesService
from scripts.KeggScraper import KeggScraper
from utils.FileReader import FileReader
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.logger import logger


@Inject(gtf_genes_service=GtfGenesService,
        file_reader=FileReader,
        whole_chromosomes_service=WholeChromosomesService,
        kegg_categories_service=KeggCategoriesService,
        kegg_subcategories_service=KeggSubcategoriesService,
        gtf_genes_kegg_categories_service=GtfGenesKeggCategoriesService,
        gtf_genes_kegg_subcategories_service=GtfGenesKeggSubcategoriesService,
        kegg_scraper=KeggScraper)
class GenesLoader:

    def __init__(self, whole_chromosomes_service=None, gtf_genes_service=None, file_reader=None,
                 kegg_categories_service=None, kegg_subcategories_service=None, gtf_genes_kegg_categories_service=None,
                 gtf_genes_kegg_subcategories_service=None, kegg_scraper=None):
        self.whole_chromosomes_service = whole_chromosomes_service
        self.gtf_genes_service = gtf_genes_service
        self.kegg_categories_service = kegg_categories_service
        self.kegg_subcategories_service = kegg_subcategories_service
        self.gtf_genes_categories_service = gtf_genes_kegg_categories_service
        self.gtf_genes_subcategories_service = gtf_genes_kegg_subcategories_service
        self.file_reader = file_reader
        self.kegg_scraper = kegg_scraper

    @DBConnection
    @TryExcept
    @Timer
    def load_genes_from_file(self, path):
        logger.info(f"Loading genes data to database from {path}")
        gene_df = self._read_and_filter_gene_data(path)
        chromosome_dfs = self.file_reader.divide_genome_df_rows_by_chromosome(gene_df)

        for chromosome_df in chromosome_dfs:
            self._process_chromosome_data(chromosome_df)

    def _read_and_filter_gene_data(self, path):
        file_df = self.file_reader.read_gtf_file(path)
        return file_df[file_df['feature'] == 'gene'].reset_index(drop=True)

    def _process_chromosome_data(self, chromosome_df):
        for _, row in chromosome_df.iterrows():
            self._process_gene_row(row)

    def _process_gene_row(self, row):
        try:
            whole_chromosome_id = self.whole_chromosomes_service.extract_id_by_refseq_accession_number(
                row['refseq_accession_number'])
        except Exception as e:
            logger.error(
                f"Failed to extract whole chromosome ID for refseq_accession_number {row['refseq_accession_number']}: {e}")
            return

        attributes = self.file_reader.divide_gtf_attributes(row['attributes'])
        record = self._build_record(row, whole_chromosome_id, attributes)
        gtf_genes_id = self.gtf_genes_service.insert(record=record)
        gene_id_gtf = attributes.get('gene_id', "")
        #self.__add_categories_and_subcategories_to_gene(gtf_genes_id, gene_id_gtf)


    def _build_record(self, row, whole_chromosome_id, attributes):

        return (
            whole_chromosome_id,
            row['source'],
            row['feature'],
            row['start_position'],
            row['end_position'],
            row['length'],
            row['score'],
            row['strand'],
            row['frame'],
            attributes.get('gene_id', ""),
            attributes.get('gene', ""),
            attributes.get('gene_biotype', "")
        )

    def __add_categories_and_subcategories_to_gene(self, gtf_genes_id: int, gene_id_gtf: str):
        logger.info(f"Starting adding KEGG categories and subcategories for gene {gene_id_gtf}")

        categories_result = self.kegg_scraper.get_categories_and_subcategories_by_gene_id_gtf(gene_id_gtf)
        if categories_result is None:
            logger.warning(f"No categories or subcategories identified for gene {gene_id_gtf}")
            return

        categories, subcategories = categories_result
        categories_ids = []
        subcategories_ids = []
        for category in categories:
            categories_ids.append(self.kegg_categories_service.extract_id_by_category_name(category))
        for subcategory in subcategories:
            subcategories_ids.append(self.kegg_subcategories_service.extract_id_by_subcategory_name(subcategory))
        for category_id in categories_ids:
            self.gtf_genes_categories_service.insert(record=(gtf_genes_id, category_id))
        for subcategory_id in subcategories_ids:
            self.gtf_genes_subcategories_service.insert(record=(gtf_genes_id, subcategory_id))
        logger.info(f"Finished adding KEGG categories and subcategories for gene {gene_id_gtf}")
