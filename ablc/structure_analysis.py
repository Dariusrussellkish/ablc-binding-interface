#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

from Bio.PDB import *
from Bio.PDB import Structure, Chain, Residue, Atom
import numpy as np
# from collections import namedtuple

# BoundingBox = namedtuple("BoundingBox", ("min_x", "max_x", "min_y", "max_y", "min_z", "max_z"))

pdb_parser = PDBParser()
pdb_list = PDBList()


def download_or_open_pdb(pdb_id: str, default_save_dir="tests/data/pdb/"):
    pdb_location = pdb_list.retrieve_pdb_file(pdb_id, pdir=default_save_dir, file_format="pdb")
    structure = pdb_parser.get_structure(pdb_id, pdb_location)
    assert structure is not None
    return structure


def substructures_to_atoms_list(substrs):
    arr = []
    for substr in substrs:
        arr += list(substr.get_atoms())
    return arr


def substructures_to_coordinates(substrs):
    def atom_to_coordinates(atom: Atom.Atom):
        return atom.get_coord()

    arr = []
    for substr in substrs:
        arr += list(map(atom_to_coordinates, substr.get_atoms()))
    return np.stack(arr)


def compute_distance_matrix(lc_coords: np.ndarray, antigen_coords: np.ndarray):
    assert lc_coords.shape[1] == 3 and lc_coords.ndim == 2
    assert antigen_coords.shape[1] == 3 and antigen_coords.ndim == 2

    distance_matrix = np.linalg.norm(
        lc_coords[:, None, :] - antigen_coords[None, :, :], axis=-1
    )

    return distance_matrix


def filter_distance_matrix(dist_mat: np.ndarray, cutoff: float):
    assert dist_mat.ndim == 2
    res = np.where(dist_mat <= cutoff)
    return list(zip(res[0], res[1]))


def get_residues_from_filter_results(lc_atoms, antigen_atoms, filter_results):
    residue_hit_set = set()
    for lc_atom_idx, ag_atom_idx in filter_results:
        lc_res = lc_atoms[lc_atom_idx].get_parent()
        ag_res = antigen_atoms[ag_atom_idx].get_parent()

        residue_hit_set.add((lc_res, ag_res))
    return residue_hit_set

# TODO: use a more efficient "collision" detection mechanism
#
# def bounding_box_from_coordinates(coords: np.ndarray):
#     min_x, max_x = np.min(coords[:, 0]), np.max(coords[:, 0])
#     min_y, max_y = np.min(coords[:, 1]), np.max(coords[:, 1])
#     min_z, max_z = np.min(coords[:, 2]), np.max(coords[:, 2])
#
#     return BoundingBox(min_x, max_x, min_y, max_y, min_z, max_z)
