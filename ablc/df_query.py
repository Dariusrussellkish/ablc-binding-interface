#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import pandas as pd
import numpy as np


def filter_for_unique_chains(df: pd.DataFrame, heavy_chain_header: str = 'Hchain',
                             light_chain_header: str = 'Lchain',
                             antigen_chain_header: str = 'antigen_chain'):
    mask = df[heavy_chain_header] != df[light_chain_header]
    mask &= df[heavy_chain_header] != df[antigen_chain_header]
    mask &= df[light_chain_header] != df[antigen_chain_header]

    return df[mask]


def filter_for_peptide_antigens(df: pd.DataFrame, antigen_chain_header: str = 'antigen_type',
                                kw: str = 'protein'):
    mask = df[antigen_chain_header] == kw
    return df[mask]


def filter_by_resolution(df: pd.DataFrame, cutoff: float, resolution_header: str = 'resolution'):
    mask = df[resolution_header] <= cutoff
    return df[mask]


def combine_dup_rows(df: pd.DataFrame):
    df = df.astype({'pdb': str, 'Hchain': str, 'Lchain': str, "model": int, "antigen_chain": str})
    df = df.groupby('pdb').agg({'Hchain': set,
                                'Lchain': set,
                                'antigen_chain': set,
                                'model': 'first'}).reset_index()

    mask = np.full((df.shape[0], 1), False)
    for i, row in enumerate(df.itertuples(index=False)):
        mask[i] = len(row.Hchain.intersection(row.Lchain)) == 0 and \
                  len(row.Hchain.intersection(row.antigen_chain)) == 0 and \
                  len(row.Lchain.intersection(row.antigen_chain)) == 0

    return df[mask]
