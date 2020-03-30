# -*- coding: utf-8 -*-

"""Utilities for INDRA."""

import itertools as itt
import logging
from dataclasses import dataclass
from operator import attrgetter
from typing import Collection, Iterable, List, Optional, TextIO, Union

from indra.assemblers.pybel import PybelAssembler
from indra.sources import indra_db_rest
from indra.statements import Evidence, Statement
from indra.tools.assemble_corpus import filter_belief, run_preassembly
from pybel import BELGraph
from pybel.canonicalize import edge_to_tuple
from pybel.constants import ANNOTATIONS, CITATION, CITATION_IDENTIFIER, EVIDENCE, RELATION, UNQUALIFIED_EDGES

__all__ = [
    'get_and_write_statements_from_agents',
    'get_and_write_statements_from_pmids',
    'get_rows_from_statement',
    'get_rows_from_statements',
    'get_graph_from_statement',
]

logger = logging.getLogger(__name__)

start_header = [
    'PMID',
    'Evidence',
    'Subject',
    'Predicate',
    'Object',
]
end_header = [
    'INDRA UUID',
    'Source Hash',
    'API',
    'Belief',
]


@dataclass
class Row:
    uuid: str
    evidence_source_hash: str
    api: str
    belief: float
    pmid: str
    evidence: str
    bel_subject: str
    bel_relation: str
    bel_object: str

    @property
    def start_tuple(self):
        return self.pmid, self.evidence, self.bel_subject, self.bel_relation, self.bel_object

    @property
    def end_tuple(self):
        return self.uuid, self.evidence_source_hash, self.api, self.belief


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
    file: Union[None, str, TextIO] = None,
    sep: Optional[str] = None,
    limit: Optional[int] = None,
    duplicates: bool = False,
    keep_only_query_pmids: bool = False,
    minimum_belief: Optional[float] = None,
    extra_columns: Optional[List[str]] = None,
) -> List[Statement]:
    """Get INDRA statements for the given agents and write the to a TSV for BEL curation.

    :param pmids: A finite iterable of PubMed identifiers
    :param file: The file to write to
    :param sep: The separator for the CSV. Defaults to a tab.
    :param limit: The optional limit of statements to write
    :param duplicates: should duplicate statements be written (with multiple evidences?)
    :param keep_only_query_pmids: If set only keeps evidences from this PMID. Warning: still might
     have multiple evidences.
    :param minimum_belief: The minimum belief score to keep
    :param extra_columns: Headers of extra columns for curation
    """
    if isinstance(pmids, str):
        pmids = [pmids]

    ids = [('pmid', pmid.strip()) for pmid in pmids]
    statements = indra_db_rest.get_statements_for_paper(ids=ids, simple_response=True)

    print_statements(
        statements,
        file=file,
        sep=sep,
        limit=limit,
        duplicates=duplicates,
        keep_only_pmids=pmids if keep_only_query_pmids else None,
        minimum_belief=minimum_belief,
        extra_columns=extra_columns,
    )

    return statements


def print_statements(
    statements: List[Statement],
    file: Union[None, str, TextIO] = None,
    sep: Optional[str] = None,
    limit: Optional[int] = None,
    duplicates: bool = False,
    keep_only_pmids: Union[None, str, Collection[str]] = None,
    sort_attrs: Iterable[str] = ('uuid', 'pmid'),
    minimum_belief: Optional[float] = None,
    extra_columns: Optional[List[str]] = None,
) -> None:
    """Write statements to a CSV for curation.

    This one is similar to the other one, but sorts by the BEL string and only keeps the first for each group.
    """
    sep = sep or '\t'
    extra_columns = extra_columns or []
    extra_columns_placeholders = [''] * len(extra_columns)

    statements = run_preassembly(statements)

    if minimum_belief is not None:
        statements = filter_belief(statements, minimum_belief)

    rows = get_rows_from_statements(statements, duplicates=duplicates, keep_only_pmids=keep_only_pmids)
    rows = sorted(rows, key=attrgetter(*sort_attrs))

    if limit is not None:
        rows = rows[:limit]

    if not rows:
        logger.warning('no rows written')
        return

    def _write(_file):
        print(*start_header, *extra_columns, *end_header, sep=sep, file=_file)
        for row in rows:
            print(*row.start_tuple, *extra_columns_placeholders, *row.end_tuple, sep=sep, file=_file)

    if isinstance(file, str):
        with open(file, 'w') as _file:
            _write(_file)
    else:
        _write(file)


def get_rows_from_statements(
    statements: Iterable[Statement],
    duplicates: bool = False,
    keep_only_pmids: Union[None, str, Collection[str]] = None,
) -> List[Row]:
    """Build and sort BEL curation rows from a list of statements using only the first evidence for each."""
    for statement in statements:
        yield from get_rows_from_statement(statement, duplicates=duplicates, keep_only_pmids=keep_only_pmids)


def get_rows_from_statement(
    statement: Statement,
    duplicates: bool = True,
    keep_only_pmids: Union[None, str, Collection[str]] = None,
) -> Iterable[Row]:
    """Convert an INDRA statement into an iterable of BEL curation rows.

    :param statement: The INDRA statement
    :param duplicates: Keep several evidences for the same INDRA statement
    :param keep_only_pmids: If set only keeps evidences from this PMID. Warning: still might
     have multiple evidences.
    """
    statement.evidence = [e for e in statement.evidence if _keep_evidence(e)]

    # Remove evidences from BioPax
    if 0 == len(statement.evidence):
        return iter([])

    if isinstance(keep_only_pmids, str):
        keep_only_pmids = {keep_only_pmids}
    if keep_only_pmids is not None:
        statement.evidence = [
            evidence
            for evidence in statement.evidence
            if evidence.pmid in keep_only_pmids
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

        bel_subject, bel_relation, bel_object = edge_to_tuple(u, v, data, use_identifiers=True)

        if any(
            substring in bel_part
            for bel_part, substring in itt.product((bel_subject, bel_relation, bel_object), SUBSTRING_BLACKLIST)
        ):
            continue

        if CITATION not in data:
            logger.debug('no citation information')
            continue

        yield Row(
            uuid=statement.uuid,
            evidence_source_hash=data[ANNOTATIONS]['source_hash'],
            belief=round(statement.belief, 2),
            pmid=data[CITATION][CITATION_IDENTIFIER],
            evidence=data[EVIDENCE],
            api=data[ANNOTATIONS]['source_api'],
            bel_subject=bel_subject,
            bel_relation=bel_relation,
            bel_object=bel_object,
        )


def get_graph_from_statement(statement: Statement) -> BELGraph:
    """Convert an INDRA statement to a BEL graph."""
    pba = PybelAssembler([statement], allow_ungrounded=True)

    try:
        graph = pba.make_model()
    except AttributeError:  # something funny happening
        logger.exception('problem making BEL graph from INDRA statements')
        return BELGraph()
    else:
        return graph
