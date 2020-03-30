# -*- coding: utf-8 -*-

"""The command line interface for BEL enrichment."""

import os
import pickle
import sys
from typing import BinaryIO, List, TextIO

import click

from pybel import BELGraph
from pybel.cli import graph_pickle_argument
from .indra_utils import get_and_write_statements_from_agents, get_and_write_statements_from_pmids
from .ranking import process_rank_genes
from .workflow import export_separate

info_cutoff_option = click.option(
    '--info-cutoff',
    type=float,
    default=1.0,
    help='Minimum inverse node degree. Lower allows more genes.',
)
belief_cutoff_option = click.option(
    '--belief-cutoff',
    type=float,
    default=0.30,
    help='Minimum belief score. Lower gets more statements.',
)
only_query_option = click.option('--only-query', is_flag=True)


@click.group()
def main():
    """BEL Enrichment."""


@main.command()
@graph_pickle_argument
@click.option('-n', '--number', type=int)
@click.option('-s', '--sep', default='\t')
def ranks(graph: BELGraph, number, sep):
    """Rank the genes in a graph."""
    gene_map = process_rank_genes(graph)
    for (namespace, name), rank in gene_map.most_common(n=number):
        click.echo(f'{rank:.2f}{sep}{namespace}{sep}{name}')


@main.command()
@graph_pickle_argument
@click.option('-d', '--directory', type=click.Path(file_okay=False, dir_okay=True), default=os.getcwd())
@info_cutoff_option
@belief_cutoff_option
def from_graph(graph: BELGraph, directory: str, info_cutoff: float, belief_cutoff: float):
    """Make a a sheet for rational enrichment of the given BEL graph."""
    export_separate(
        graph=graph,
        directory=directory,
        minimum_information_density=info_cutoff,
        minimum_belief=belief_cutoff,
    )


output_option = click.option('--output', type=click.File('w'), default=sys.stdout, help='output file')
pickle_output_option = click.option('--pickle-file', type=click.File('wb'), help='output file')
no_duplicates_option = click.option('--no-duplicates', is_flag=True)


@main.command()
@click.option('-a', '--agents', multiple=True)
@output_option
@pickle_output_option
@belief_cutoff_option
@no_duplicates_option
def from_agents(agents: List[str], output: TextIO, pickle_file: BinaryIO, belief_cutoff: float, no_duplicates: bool):
    """Make a sheet for the given agents."""
    statements = get_and_write_statements_from_agents(
        agents=agents,
        file=output,
        duplicates=(not no_duplicates),
        minimum_belief=belief_cutoff,
    )

    if pickle_file:
        pickle.dump(statements, pickle_file)


@main.command()
@click.argument('pmids', nargs=-1)
@output_option
@pickle_output_option
@belief_cutoff_option
@no_duplicates_option
@only_query_option
def from_pmids(
    pmids: List[str],
    output: TextIO,
    pickle_file: BinaryIO,
    belief_cutoff: float,
    no_duplicates: bool,
    only_query: bool,
):
    """Make a sheet for the given PMIDs."""
    statements = get_and_write_statements_from_pmids(
        pmids=pmids,
        file=output,
        duplicates=(not no_duplicates),
        keep_only_query_pmids=only_query,
        minimum_belief=belief_cutoff,
    )

    if pickle_file:
        pickle.dump(statements, pickle_file)


@main.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin, help='a text file with one PMID per line')
@output_option
@pickle_output_option
@belief_cutoff_option
@no_duplicates_option
@only_query_option
def from_pmid_file(
    pmids: TextIO,
    output: TextIO,
    pickle_file: BinaryIO,
    belief_cutoff: float,
    no_duplicates: bool,
    only_query: bool,
):
    """Make a sheet for the PMIDs in the given file."""
    statements = get_and_write_statements_from_pmids(
        pmids=pmids,
        file=output,
        duplicates=(not no_duplicates),
        minimum_belief=belief_cutoff,
        keep_only_query_pmids=only_query,
    )

    if pickle_file:
        pickle.dump(statements, pickle_file)


if __name__ == '__main__':
    main()
