import requests
import re
from bs4 import BeautifulSoup

from .kegg_pathway import KEGG_PATHWAY_DICT

from utils.logger import logger


class KeggScraper():
    def _get_org_code_and_gene_id(self, gene_id: str):
        url = (f"https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=genes"
               f"&keywords={gene_id}&page=1")
        response = requests.get(url)

        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            a_tags = soup.find_all('a')
            matching_links = []
            for a in a_tags:
                if a.string and gene_id in a.string:
                    matching_links.append(a)
            if len(matching_links) == 0:
                return None
            return matching_links[0].get_text(strip=True)
        else:
            raise Exception(f"Failed to retrieve the page. Status code: {response.status_code}")

    def _get_ko_url(self, org_code_and_gene_id: str):
        if org_code_and_gene_id is None:
            return None
        url = f"https://www.genome.jp/entry/{org_code_and_gene_id}"
        response = requests.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            rows = soup.find_all('tr')
            for row in rows:
                ko_span = row.find('span')
                if ko_span and ko_span.string:
                    if 'KO' in ko_span.string:
                        a_tag = row.find('a')
                        if a_tag:
                            return a_tag.get_text(strip=True)
            logger.warning(f"No KO entry found for gene - {org_code_and_gene_id}")
            return None
        else:
            raise Exception(f"Failed to retrieve the page. Status code: {response.status_code}")

    def _get_orthology_text(self, ko: str):
        if ko is None:
            return None
        url = f"https://www.genome.jp/entry/{ko}"
        response = requests.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            rows = soup.find_all('tr')
            for row in rows:
                th = row.find('th')
                if th and th.get_text(strip=True) == "Brite":
                    cel_div = row.find('div', class_='cel')
                    if cel_div:
                        raw_text = cel_div.get_text(strip=False)
                        cleaned_text = re.sub(r'\s+', ' ', raw_text)
                        formatted_text = re.sub(r'(\D)(K\d+|\d+)', r'\1 \2', cleaned_text)
                        return formatted_text

            logger.warning("No 'Brite' entry found with the associated 'cel' div.")
            return None
        else:
            raise Exception(f"Failed to retrieve the page. Status code: {response.status_code}")

    def _get_categories_with_complete_dict(self, categories, orthology_text):
        if orthology_text is None:
            return None
        identified_categories = []
        identified_subcategories = []
        identified_subsubcategories = []
        for category, subcategories_dict in categories.items():
            if category.lower() in orthology_text:
                identified_categories.append(category)
            for subcategory, subsubcategories_list in subcategories_dict.items():
                if subcategory.lower() in orthology_text:
                    identified_subcategories.append(subcategory)
                    identified_categories.append(category)
                for subsubcategory in subsubcategories_list:
                    if subsubcategory.lower() in orthology_text:
                        identified_subsubcategories.append(subsubcategory)
                        identified_subcategories.append(subcategory)
                        identified_categories.append(category)
        return list(set(identified_categories)), list(set(identified_subcategories))

    def get_categories_and_subcategories_by_gene_id_gtf(self, gene_id_gtf: str) -> tuple[list[str], list[str]] | None:
        text = self._get_orthology_text(self._get_ko_url(self._get_org_code_and_gene_id(gene_id_gtf)))
        text = text.lower() if text else text
        categories_result = self._get_categories_with_complete_dict(KEGG_PATHWAY_DICT, text)
        return categories_result
