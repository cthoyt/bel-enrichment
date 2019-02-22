# -*- coding: utf-8 -*-

"""Load a BEL graph from curation sheets."""

import logging
import os
from collections import defaultdict
from typing import Dict, Iterable, Mapping, Optional, Tuple

import pandas as pd
import pyparsing
from tqdm import tqdm

from pybel.constants import CITATION_REFERENCE, CITATION_TYPE, CITATION_TYPE_PUBMED
from pybel.parser import BELParser
from pybel.parser.exc import BELParserWarning, BELSyntaxError

logger = logging.getLogger(__name__)

NOT_CURATED = 'Not curated'
ERROR = 'Error'
CORRECT = 'Correct'
ERROR_BUT_ALSO_OTHER_STATEMENT = 'Error but other statement was identified'
MODIFIED_BY_CURATOR = 'Modified by curator'


def _check_curation_template_columns(df: pd.DataFrame, path: str) -> bool:
    """Check the columns in a curation dataframe."""
    if 'Curator' not in df.columns:
        logger.warning(f'{path} is missing the "Curator" column')
        return False

    if 'Checked' not in df.columns:
        logger.warning(f'{path} is missing the "Checked" column')
        return False

    if 'Correct' not in df.columns:
        logger.warning(f'{path} is missing the "Correct" column')
        return False

    if 'Changed' not in df.columns:
        logger.warning(f'{path} is missing the "Changed" column')
        return False

    return True


def process_row(bel_parser: BELParser, row: Dict, line_number: int) -> None:
    """Process a row."""
    if not row['Checked']:  # don't use unchecked material
        return

    if not (row['Correct'] or row['Changed']):  # if it's neither correct nor changed, then it's fucked
        return

    reference = str(
        row['Citation Reference']
        if 'Citation Reference' in row else
        row['PMID']
    )

    if not reference:
        raise Exception('missing reference')

    bel_parser.control_parser.citation = {
        CITATION_TYPE: CITATION_TYPE_PUBMED,
        CITATION_REFERENCE: reference,
    }
    # Set the evidence
    bel_parser.control_parser.evidence = row['Evidence']
    # TODO set annotations if they exist

    annotations = {
        'Curator': row['Curator'],
        'Confidence': 'Medium',  # needs re-curation
    }

    if 'INDRA UUID' in row:
        annotations['INDRA_UUID'] = row['INDRA UUID']

    if 'Belief' in row:
        annotations['INDRA_Belief'] = row['Belief']

    if 'API' in row:
        annotations['INDRA_API'] = row['API']

    # Set annotations
    bel_parser.control_parser.annotations.update(annotations)

    sub = row['Subject']
    obj = row['Object']

    # Build a BEL statement and parse it
    bel = f"{sub} {row['Predicate']} {obj}"

    # Cast line number from numpy.int64 to integer since JSON cannot handle this class
    line_number = int(line_number)

    try:
        bel_parser.parseString(bel, line_number=line_number)
    except BELParserWarning as exc:
        bel_parser.graph.add_warning(exc)
    except pyparsing.ParseException as exc:
        bel_parser.graph.add_warning(BELSyntaxError(line_number=line_number, line=bel, position=exc.loc))


def generate_error_types(path: str) -> Tuple[Mapping[str, int], str]:
    """Generate report about the types of errors INDRA made.

    :param path: path to the excel file
    :return: summary of the curation
    """
    df = pd.read_excel(path)

    error_types = defaultdict(int)
    curator = None

    for line, row in df.iterrows():
        if line == 0:
            curator = row.get('Curator')

        error_type = row.get('Error Type')

        if pd.isnull(error_type):
            continue

        # Multiple errors are listed separated by a comma
        for error in str(error_type).split(','):
            # Lower case errors and remove spacing
            error_types[error.lower().strip()] += 1

    return error_types, curator


def generate_curation_report(path: str) -> dict:
    """Generate report about curated/non-curated statements in a given curation template.

    :param path: path to the excel file
    :return: summary of the curation
    """
    df = pd.read_excel(path)

    # Check columns in dataframe exist
    if not _check_curation_template_columns(df, path):
        raise ValueError(f'{path} has a problem with the header')

    curation_results = defaultdict(int)

    for line, row in df.iterrows():
        evidence = row.get('Evidence')
        if evidence == 'No evidence text.':
            logger.debug('No evidence text. Skipping...')
            continue

        checked = row.get('Checked')
        correct = row.get('Correct')
        changed = row.get('Changed')
        # Transform real values ('x' and'NaN') to Trues and Falses
        checked = pd.notnull(checked)
        correct = pd.notnull(correct)
        changed = pd.notnull(changed)

        # The statement has not been curated (all 3 columns are empty)
        if not any([checked, correct, changed]):
            curation_results[NOT_CURATED] += 1

        # Only checked is marked
        elif checked and not any([correct, changed]):
            curation_results[ERROR] += 1

        # Correct statements by Indra
        elif correct and not changed:
            curation_results[CORRECT] += 1

        # Statement has been modified by the curator and WAS the original one
        elif checked and changed:
            curation_results[MODIFIED_BY_CURATOR] += 1

        elif changed and correct:
            logger.warning(f'Conflict in row {line}')

        # Statement has been modified by the curator but WAS NOT the original one
        elif changed:
            curation_results[ERROR_BUT_ALSO_OTHER_STATEMENT] += 1

        curation_results['Total'] += 1

    return dict(curation_results)


def generate_curation_summary(input_directory: str,
                              output_directory: str,
                              use_tqdm: bool = True,
                              sheet_suffix: Optional[str] = None,
                              ) -> None:
    """Generate a summary of the curation results on excel."""
    summary_excel_rows = {}
    error_excel_rows = {}

    paths = iterate_sheets_paths(input_directory, suffix=sheet_suffix)
    if use_tqdm:
        paths = tqdm(list(paths), desc='Generating curation report')

    for path in paths:
        gene_symbol = path.split('/')[-2]

        # Subfolder name (Gene Symbol) -> dictionary results
        summary_excel_rows[gene_symbol] = generate_curation_report(path)

        error_types, _ = generate_error_types(path)

        error_excel_rows[gene_symbol] = error_types

    # Export Summary Report
    df_summary = pd.DataFrame.from_dict(summary_excel_rows, orient='index')
    df_summary = df_summary.fillna(0).astype(int)
    # Rearrange columns
    df_summary = df_summary[[CORRECT, ERROR, ERROR_BUT_ALSO_OTHER_STATEMENT, MODIFIED_BY_CURATOR, NOT_CURATED, 'Total']]
    df_summary.to_csv(os.path.join(output_directory, 'curation_summary.csv'))

    # Export Error Types Report
    df_error = pd.DataFrame.from_dict(error_excel_rows, orient='index')
    df_error = df_error.fillna(0).astype(int)
    df_error.to_csv(os.path.join(output_directory, 'error_types.csv'))


def iterate_sheets_paths(directory: str, suffix: Optional[str] = None) -> Iterable[str]:
    """List the excel curation sheets."""
    if suffix is None:
        suffix = '_curated.xlsx'
    for path in os.listdir(directory):
        subpath = os.path.join(directory, path)
        if os.path.isdir(subpath):
            yield from iterate_sheets_paths(subpath, suffix=suffix)
        elif subpath.endswith(suffix):
            yield subpath
