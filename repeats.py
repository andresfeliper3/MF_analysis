from src.Biocode.services.RMRepeatsWholeChromosomesService import RMRepeatsWholeChromosomesService
from src.Biocode.services.RepeatsService import RepeatsService
from src.Biocode.services.WholeChromosomesService import WholeChromosomesService

from load import loader
from utils.decorators import Timer, DBConnection
from utils.logger import logger

import pandas as pd

def _import_file_out(path):
  data = []

  with open(path, 'r') as file:
      next(file)  # Skip 3 lines
      next(file)
      next(file)
      for line in file:
          # Split the line by whitespace(s)
          columns = line.strip().split()
          data.append(columns)

  df = pd.DataFrame(data, columns=[
    'sw_score', 'percentage_divergence', 'percentage_deletions', 'percentage_insertions', 'sequence', 'query_begin',
    'query_end', 'query_left', 'strand', 'name', 'class_family', 'repeat_begin', 'repeat_end', 'repeat_left', 'ID', 'add'])

  df['query_end'] = df['query_end'].astype('int')
  df['query_begin'] = df['query_begin'].astype('int')
  df["repeat_length"] = df["query_end"] - df["query_begin"] + 1
  return df

@DBConnection
@Timer
def load_RM_repeats_from_file(path):
    logger.info(f"Loading repeats results from RepeatMasker from {path}")
    df = _import_file_out(path)
    repeats_service = RepeatsService()
    rmRepeatsWholeChromosomesService= RMRepeatsWholeChromosomesService()
    wholeChromosomesService = WholeChromosomesService()

    for _, row in df.iterrows():
        repeat_id = repeats_service.insert(record=(row['name'], row['class_family'], 'RM'))
        whole_chromosome_id = wholeChromosomesService.extract_id_by_refseq_accession_number(row['sequence'])
        record = (
            repeat_id,
            whole_chromosome_id,
            row['sw_score'],
            row['percentage_divergence'],
            row['percentage_deletions'],
            row['percentage_insertions'],
            row['query_begin'],
            row['query_end'],
            row['repeat_length'],
            row['query_left'],
            row['strand'],
            row['repeat_begin'],
            row['repeat_end'],
            row['repeat_left']
        )
        rmRepeatsWholeChromosomesService.insert(record=record)

