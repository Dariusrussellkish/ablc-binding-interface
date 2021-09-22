#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

from Bio.PDB import *
from Bio.PDB import Atom, Chain, Residue
import numpy as np
import typing as t
from .StructSeqRecord import BindingInterface


pdb_parser = PDBParser(QUIET=True)
pdb_list = PDBList(verbose=False)


def substructures_to_residues_list(substrs: t.Collection) -> t.List[Residue.Residue]:
    arr = []
    for substr in substrs:
        arr += list(substr.get_residues())
    return arr


def substructures_to_atoms_list(substrs: t.Collection) -> t.List[Atom.Atom]:
    arr = []
    for substr in substrs:
        arr += list(substr.get_atoms())
    return arr


def substructures_to_coordinates(substrs: t.Collection) -> np.ndarray:
    def atom_to_coordinates(atom: Atom.Atom):
        return atom.get_coord()

    arr = []
    for substr in substrs:
        arr += list(map(atom_to_coordinates, substr.get_atoms()))
    return np.vstack(arr)


def compute_distance_matrix(lc_coords: np.ndarray,
                            antigen_coords: np.ndarray) -> np.ndarray:
    assert lc_coords.shape[1] == 3 and lc_coords.ndim == 2
    assert antigen_coords.shape[1] == 3 and antigen_coords.ndim == 2

    distance_matrix = np.linalg.norm(
        lc_coords[:, None, :] - antigen_coords[None, :, :], axis=-1
    )

    return distance_matrix


def filter_distance_matrix(dist_mat: np.ndarray,
                           cutoff: float) -> t.List[t.Tuple[int, int]]:
    assert dist_mat.ndim == 2
    res = np.where(dist_mat <= cutoff)
    return list(zip(res[0], res[1]))


def get_residues_from_filter_results(lc_atoms: t.Sequence[Atom.Atom],
                                     antigen_atoms: t.Sequence[Atom.Atom],
                                     filter_results) -> t.Set[t.Tuple[Residue.Residue, Residue.Residue]]:
    residue_hit_set = set()
    for lc_atom_idx, ag_atom_idx in filter_results:
        lc_res = lc_atoms[lc_atom_idx].get_parent()
        ag_res = antigen_atoms[ag_atom_idx].get_parent()

        if lc_res == ag_res:
            continue

        residue_hit_set.add((lc_res, ag_res))
    return residue_hit_set


def get_binding_partners(a: t.Set[Chain.Chain],
                         b: t.Set[Chain.Chain],
                         combine=False,
                         cutoff=6.0) -> BindingInterface:
    """
    Determine the binding interface residues of two sets of chains

    Note: it will ignore identical residue pairs on an interface caused by chains in the intersection of a and b.
    Normally, the intersection should be empty

    :param a: Set of chain names for a
    :param b: Set of chain names for b
    :param combine: Combine a and b to do an all-to-all comparison
    :param cutoff: cutoff distance to consider for contacts
    :return: BindingInterface
    """

    a_chains = {c.get_id() for c in a}
    b_chains = {c.get_id() for c in b}

    a_atoms = substructures_to_atoms_list(a)
    b_atoms = substructures_to_atoms_list(b)

    a_coords = substructures_to_coordinates(a)
    b_coords = substructures_to_coordinates(b)

    if combine:
        atoms = a_atoms + b_atoms
        coords = np.vstack([a_coords, b_coords])

        dist = compute_distance_matrix(coords, coords)
        indices = filter_distance_matrix(dist, cutoff)
        return BindingInterface(
            a_chains,
            b_chains,
            get_residues_from_filter_results(atoms, atoms, indices)
            )
    else:
        dist = compute_distance_matrix(a_coords, b_coords)
        indices = filter_distance_matrix(dist, cutoff)
        return BindingInterface(
            a_chains,
            b_chains,
            get_residues_from_filter_results(a_atoms, b_atoms, indices)
        )
