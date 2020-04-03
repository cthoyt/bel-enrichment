# -*- coding: utf-8 -*-

"""The command line interface for BEL enrichment."""

import json
import os
import sys
from typing import List, TextIO

import click

from pybel import BELGraph
from pybel.cli import graph_pickle_argument

from indra.statements import stmts_to_json
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
statement_json_file_option = click.option('--statement-file', type=click.File('w'), help='output statements JSON file')
no_duplicates_option = click.option('--no-duplicates', is_flag=True)
no_ungrounded_option = click.option('--no-ungrounded', is_flag=True)


@main.command()
@click.option('-a', '--agents', multiple=True)
@output_option
@statement_json_file_option
@belief_cutoff_option
@no_duplicates_option
@no_ungrounded_option
def from_agents(
    agents: List[str],
    output: TextIO,
    statement_file: TextIO,
    belief_cutoff: float,
    no_duplicates: bool,
    no_ungrounded: bool,
):
    """Make a sheet for the given agents."""
    statements = get_and_write_statements_from_agents(
        agents=agents,
        file=output,
        allow_duplicates=(not no_duplicates),
        allow_ungrounded=(not no_ungrounded),
        minimum_belief=belief_cutoff,
    )

    if statement_file:
        json.dump(stmts_to_json(statements), statement_file, indent=2)


@main.command()
@click.argument('pmids', nargs=-1)
@output_option
@statement_json_file_option
@belief_cutoff_option
@no_duplicates_option
@only_query_option
def from_pmids(
    pmids: List[str],
    output: TextIO,
    statement_file: TextIO,
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

    if statement_file:
        json.dump(stmts_to_json(statements), statement_file, indent=2)


@main.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin, help='a text file with one PMID per line')
@output_option
@statement_json_file_option
@belief_cutoff_option
@no_duplicates_option
@only_query_option
def from_pmid_file(
    pmids: TextIO,
    output: TextIO,
    statement_file: TextIO,
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

    if statement_file:
        json.dump(stmts_to_json(statements), statement_file, indent=2)


if __name__ == '__main__':
    main()
