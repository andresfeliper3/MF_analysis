# This script only needs to be executed ONCE before starting loading gtf_genes

from scripts.kegg_pathway import KEGG_PATHWAY_DICT
from src.Biocode.services.KeggCategoriesService import KeggCategoriesService
from src.Biocode.services.KeggSubcategoriesService import KeggSubcategoriesService
from utils.decorators import Timer, DBConnection, TryExcept


@DBConnection
@TryExcept
@Timer
def load_categories_tables():
    kegg_categories_service = KeggCategoriesService()
    kegg_subcategories_service = KeggSubcategoriesService()

    for category, subcategories_dict in KEGG_PATHWAY_DICT.items():
        record = (category, )
        kegg_category_id = kegg_categories_service.insert(record)
        for subcategory, subsubcategories_list in subcategories_dict.items():
            kegg_subcategory_id = kegg_subcategories_service.insert(record=(subcategory, kegg_category_id))


load_categories_tables()
