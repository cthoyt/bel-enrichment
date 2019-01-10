# -*- coding: utf-8 -*-

"""A wrapper for functions on a BEL sheets repository."""

import logging
import os
from dataclasses import dataclass
from typing import Any, Mapping, Optional

import pandas as pd
from tqdm import tqdm

from pybel import BELGraph, from_pickle, to_json_path, to_pickle
from pybel.parser import BELParser
from .sheets import _check_curation_template_columns, generate_curation_summary, get_sheets_paths, process_row

__all__ = [
    'BELSheetsRepository',
]

logger = logging.getLogger(__name__)


@dataclass
class BELSheetsRepository:
    """A repository of BEL curation sheets."""

    directory: str
    output_directory: str
    metadata: Optional[Mapping[str, str]] = None

    pickle_name = 'sheets.bel.pickle'
    json_name = 'sheets.bel.json'

    @property
    def _cache_pickle_path(self):
        return os.path.join(self.output_directory, self.pickle_name)

    def get_graph(self,
                  use_cached: bool = False,
                  use_tqdm: bool = True,
                  tqdm_kwargs: Optional[Mapping[str, Any]] = None,
                  ) -> BELGraph:
        """Get the BEL graph from all sheets in this repository.

        .. warning:: This BEL graph isn't pre-filled with namespace and annotation URLs.
        """
        if use_cached and self._cache_pickle_path is not None and os.path.exists(self._cache_pickle_path):
            return from_pickle(self._cache_pickle_path)

        graph = BELGraph(**(self.metadata or {}))
        logger.info('streamlining parser')
        bel_parser = BELParser(graph)

        paths = get_sheets_paths(self.directory)

        if use_tqdm:
            _tqdm_kwargs = dict(desc='Sheets')
            if tqdm_kwargs:
                _tqdm_kwargs.update(tqdm_kwargs)
            paths = tqdm(list(paths), **_tqdm_kwargs)

        for path in paths:
            df = pd.read_excel(path)

            # Check columns in DataFrame exist
            if not _check_curation_template_columns(df, path):
                continue

            graph.path = path

            line_rows = df.iterrows()
            # if use_tqdm:
            #    line_rows = tqdm(line_rows, total=len(df.index), leave=False, desc=path.split('/')[-1].split('_')[0])
            for line_number, row in line_rows:
                process_row(graph, bel_parser=bel_parser, row=row, line_number=line_number)

        to_pickle(graph, self._cache_pickle_path)
        to_json_path(graph, os.path.join(self.output_directory, self.json_name))
        return graph

    def generate_curation_summary(self):
        """Generate a curation summary."""
        return generate_curation_summary(
            input_directory=self.directory,
            output_directory=self.output_directory,
        )


@dataclass
class BELEnrichment:
    """A class containing all of the information necessary for the enrichment workflow."""

    prior: BELGraph
    directory: str
