# -*- coding: utf-8 -*-

"""Utilities for INDRA."""

import os
import pickle
from typing import List, Optional, TextIO

from indra.statements import Statement
from pybel import BELGraph
from .indra_utils import get_and_write_statements_from_agents
from .ranking import process_rank_genes

__all__ = [
    'export_separate',
    'export_single',
]


def export_separate(
    graph: BELGraph,
    directory: str,
    minimum_information_density: float = 1.0,
    minimum_belief: float = 0.3,
    sep: str = '\t',
    limit: Optional[int] = None,
    duplicates: bool = False,
):
    """Get genes from the graph and export in separate folders."""
    gene_symbols = get_gene_symbols(
        graph=graph,
        cutoff=minimum_information_density,
    )

    for gene_symbol in gene_symbols:
        gene_directory = os.path.join(directory, gene_symbol)
        os.makedirs(gene_directory, exist_ok=True)
        tsv_path = os.path.join(gene_directory, f'{gene_symbol}.bel.tsv')
        pickle_path = os.path.join(gene_directory, f'{gene_symbol}_statements.pkl')

        if os.path.exists(tsv_path):
            continue  # already downloaded

        with open(tsv_path, 'w') as csv_file:
            statements = get_and_write_statements_from_agents(
                agents=gene_symbol,
                file=csv_file,
                sep=sep,
                limit=limit,
                allow_duplicates=duplicates,
                minimum_belief=minimum_belief,
            )
        with open(pickle_path, 'wb') as pkl_file:
            pickle.dump(statements, pkl_file)


def export_single(
    graph: BELGraph,
    cutoff: float = 1.0,
    file: Optional[TextIO] = None,
    sep: str = '\t',
    limit: Optional[int] = None,
    duplicates: bool = False
) -> List[Statement]:
    """Get genes from the graph and export as one file."""
    gene_symbols = get_gene_symbols(graph=graph, cutoff=cutoff)

    return get_and_write_statements_from_agents(
        agents=gene_symbols,
        file=file,
        sep=sep,
        limit=limit,
        allow_duplicates=duplicates,
    )


def get_gene_symbols(graph: BELGraph, cutoff: float = 1.0):
    """Get HGNC gene symbols having above a given cutoff."""
    gene_map = process_rank_genes(graph)

    return [
        name
        for (namespace, name), rank in gene_map.most_common()
        if namespace.lower() == 'hgnc' or cutoff < rank
    ]
