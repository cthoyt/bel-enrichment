# -*- coding: utf-8 -*-

"""Utilities for INDRA."""

import logging
from collections import namedtuple
from operator import attrgetter
from typing import Iterable, List, Optional, TextIO, Tuple, Union

from indra.assemblers.pybel import PybelAssembler
from indra.sources import indra_db_rest
from indra.statements import Evidence, Statement
from indra.tools.assemble_corpus import filter_belief, run_preassembly
from pybel import BELGraph
from pybel.constants import ANNOTATIONS, CITATION, CITATION_IDENTIFIER, EVIDENCE, RELATION, UNQUALIFIED_EDGES

__all__ = [
    'get_and_write_statements_from_agents',
    'get_and_write_statements_from_pmids',
    'get_rows_from_statement',
    'get_rows_from_statements',
    'get_graph_from_statement',
]

logging.getLogger('pybel_assembler').setLevel(logging.ERROR)
logging.getLogger('assemble_corpus').setLevel(logging.ERROR)
logging.getLogger('preassembler').setLevel(logging.ERROR)

header = [
    'INDRA UUID',
    'Belief',
    'PMID',
    'Evidence',
    'API',
    'Subject',
    'Predicate',
    'Object',
]

Row = namedtuple('Row', ['uuid', 'evidence_source_hash', 'belief', 'pmid', 'evidence', 'api', 'bel'])

NO_EVIDENCE_TEXT = 'No evidence text.'
MODIFIED_ASSERTION = 'Modified assertion'
TEXT_BLACKLIST = {NO_EVIDENCE_TEXT, MODIFIED_ASSERTION}

SOURCE_BLACKLIST = {'bel', 'signor'}
SUBSTRING_BLACKLIST = {'CHEBI', 'PUBCHEM', 'act(bp('}


def get_and_write_statements_from_agents(
    agents: Union[str, List[str]],
    file: Optional[TextIO] = None,
    sep: Optional[str] = None,
    limit: Optional[int] = None,
    duplicates: bool = False,
    minimum_belief: Optional[float] = None,
) -> List[Statement]:
    """Get INDRA statements for the given agents and write the to a TSV for BEL curation.

    :param agents: A list of agents (HGNC gene symbols)
    :param file: The file to write to
    :param sep: The separator for the CSV. Defaults to a tab.
    :param limit: The optional limit of statements to write
    :param duplicates: should duplicate statements be written (with multiple evidences?)
    :param minimum_belief: The minimum belief score to keep
    """
    if isinstance(agents, str):
        agents = [agents]

    processor = indra_db_rest.get_statements(agents=agents)
    statements = processor.statements

    print_statements(
        statements,
        file=file,
        sep=sep,
        limit=limit,
        duplicates=duplicates,
        minimum_belief=minimum_belief,
    )

    return statements


def get_and_write_statements_from_pmids(
    pmids: Union[str, Iterable[str]],
    file: Optional[TextIO] = None,
    sep: Optional[str] = None,
    limit: Optional[int] = None,
    duplicates: bool = False,
    keep_only_pmid: Optional[str] = None,
    minimum_belief: Optional[float] = None,
) -> List[Statement]:
    """Get INDRA statements for the given agents and write the to a TSV for BEL curation.

    :param pmids: A finite iterable of PubMed identifiers
    :param file: The file to write to
    :param sep: The separator for the CSV. Defaults to a tab.
    :param limit: The optional limit of statements to write
    :param duplicates: should duplicate statements be written (with multiple evidences?)
    :param keep_only_pmid: If set only keeps evidences from this PMID. Warning: still might
     have multiple evidences.
    :param minimum_belief: The minimum belief score to keep
    """
    if isinstance(pmids, str):
        ids = [('pmid', pmids.strip())]
    else:
        ids = [('pmid', pmid.strip()) for pmid in pmids]

    statements = indra_db_rest.get_statements_for_paper(ids=ids, simple_response=True)

    print_statements(
        statements,
        file=file,
        sep=sep,
        limit=limit,
        duplicates=duplicates,
        keep_only_pmid=keep_only_pmid,
        minimum_belief=minimum_belief,
    )

    return statements


def print_statements(
    statements: List[Statement],
    file: Optional[TextIO] = None,
    sep: Optional[str] = None,
    limit: Optional[int] = None,
    duplicates: bool = False,
    keep_only_pmid: Optional[str] = None,
    sort_attrs: Iterable[str] = ('uuid', 'pmid'),
    minimum_belief: Optional[float] = None,
) -> None:
    """Write statements to a CSV for curation.

    This one is similar to the other one, but sorts by the BEL string and only keeps the first for each group.
    """
    sep = sep or '\t'

    print(*header, sep=sep, file=file)

    statements = run_preassembly(statements)

    if minimum_belief is not None:
        statements = filter_belief(statements, minimum_belief)

    rows = get_rows_from_statements(statements, duplicates=duplicates, keep_only_pmid=keep_only_pmid)
    rows = sorted(rows, key=attrgetter(*sort_attrs))

    if limit is not None:
        rows = rows[:limit]

    for row in rows:
        print(*row, sep=sep, file=file)


def get_rows_from_statements(
    statements: Iterable[Statement],
    duplicates: bool = False,
    keep_only_pmid: Optional[str] = None,
) -> List[Row]:
    """Build and sort BEL curation rows from a list of statements using only the first evidence for each."""
    for statement in statements:
        yield from get_rows_from_statement(statement, duplicates=duplicates, keep_only_pmid=keep_only_pmid)


def get_rows_from_statement(
    statement: Statement,
    duplicates: bool = True,
    keep_only_pmid: Optional[str] = None,
) -> Iterable[Row]:
    """Convert an INDRA statement into an iterable of BEL curation rows.

    :param statement: The INDRA statement
    :param duplicates: Keep several evidences for the same INDRA statement
    :param keep_only_pmid: If set only keeps evidences from this PMID. Warning: still might
     have multiple evidences.
    """
    statement.evidence = [e for e in statement.evidence if _keep_evidence(e)]

    # Remove evidences from BioPax
    if 0 == len(statement.evidence):
        return iter([])

    if keep_only_pmid is not None:
        statement.evidence = [
            evidence
            for evidence in statement.evidence
            if evidence.pmid == keep_only_pmid
        ]
        # Might also be a case where several evidences from
        # same document exist, but we really only want one.
    if not duplicates:
        # Remove all but the first remaining evidence for the statement
        # unused_evidences = statement.evidence[1:]
        del statement.evidence[1:]

    yield from _get_rows_from_statement(statement)


def _keep_evidence(evidence: Evidence):
    return (
        evidence.pmid and
        evidence.text and
        evidence.text not in TEXT_BLACKLIST and
        evidence.source_api and
        evidence.source_api not in SOURCE_BLACKLIST
    )


def _get_rows_from_statement(statement: Statement) -> Iterable[Row]:
    """Build a BEL graph from the given INDRA statement and iterate over rows of all possible BEL edges."""
    graph = get_graph_from_statement(statement)

    for u, v, data in graph.edges(data=True):
        if data[RELATION] in UNQUALIFIED_EDGES:
            continue

        bel = graph.edge_to_bel(u, v, data, sep='\t', use_identifiers=True)

        if any(s in bel for s in SUBSTRING_BLACKLIST):
            continue

        yield Row(
            uuid=statement.uuid,
            evidence_source_hash= data[ANNOTATIONS]['source_hash'],
            belief=round(statement.belief, 2),
            pmid=data[CITATION][CITATION_IDENTIFIER],
            evidence=data[EVIDENCE],
            api=data[ANNOTATIONS]['source_api'],
            bel=bel,
        )


def get_graph_from_statement(statement: Statement) -> BELGraph:
    """Convert an INDRA statement to a BEL graph."""
    pba = PybelAssembler([statement])

    try:
        graph = pba.make_model()
    except AttributeError:  # something funny happening
        return BELGraph()
    else:
        return graph
