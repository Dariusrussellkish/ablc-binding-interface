#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import typing as t
from collections import defaultdict

from Bio.PDB.Residue import Residue


def dict_from_resid_pairs(pairs: t.Collection[t.Tuple[Residue, Residue]],
                          second: bool = False) -> t.Dict[Residue, t.Set[Residue]]:
    resids = defaultdict(set)
    for a, b in pairs:
        if second:
            resids[b].add(a)
        else:
            resids[a].add(b)
    return dict(resids)


class BindingInterface:
    def __init__(self,
                 chains_a: t.Set[str],
                 chains_b: t.Set[str],
                 interface_pairs: t.Collection[t.Tuple[Residue, Residue]]
                 ):
        self.chains_a = chains_a
        self.chains_b = chains_b
        self.interface_pairs = interface_pairs
        self.a_dict = dict_from_resid_pairs(interface_pairs)
        self.b_dict = dict_from_resid_pairs(interface_pairs, second=True)
