from src.Biocode.services.GtfGenesService import GtfGenesService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService
from utils.FileReader import FileReader
from utils.decorators import Timer, DBConnection, TryExcept, Inject
from utils.logger import logger


@Inject(gtf_genes_service = GtfGenesService,
        file_reader = FileReader,
        whole_chromosomes_service = WholeChromosomesService)
class GenesLoader:

    def __init__(self, whole_chromosomes_service=None, gtf_genes_service=None, file_reader=None):
        self.whole_chromosomes_service = whole_chromosomes_service
        self.gtf_genes_service = gtf_genes_service
        self.file_reader = file_reader

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
        self.gtf_genes_service.insert(record=record)

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
