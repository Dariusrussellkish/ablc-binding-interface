#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import os
import pathlib
import pickle

from ablc.structure_analysis import *
from ablc.df_query import *

from Bio.PDB.StructureBuilder import PDBConstructionWarning
import warnings
warnings.simplefilter('ignore', category=PDBConstructionWarning)


current_dir = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
data_path = current_dir / "data"

df = pd.read_csv(data_path / "sabdab_nr_h_l_bound_summary.tsv", sep='\t')

df = filter_for_unique_chains(df)
df = filter_for_peptide_antigens(df)
df = combine_dup_rows(df)

data = dict()
rows = df.shape[0]
for i, row in enumerate(df.itertuples(index=False)):

    print(f"Processing: {row.pdb} | {i+1} / {rows}")
    try:
        pdb = download_or_open_pdb(row.pdb, default_save_dir=data_path / "pdb")
        model = pdb[row.model]
        l_chain_set = {model[lchain]for lchain in row.Lchain}
        ag_chain_set = {model[ag_chain]for ag_chain in row.antigen_chain}

        l_chain_atoms = substructures_to_atoms_list(l_chain_set)
        ag_chain_atoms = substructures_to_atoms_list(ag_chain_set)

        l_chain_coords = substructures_to_coordinates(l_chain_set)
        ag_chain_coords = substructures_to_coordinates(ag_chain_set)
        dist_mat = compute_distance_matrix(l_chain_coords, ag_chain_coords)

        indices = filter_distance_matrix(dist_mat, cutoff=4.0)
        binding_interface = get_residues_from_filter_results(l_chain_atoms, ag_chain_atoms, indices)

        res = (pdb, row.model, row.Lchain, row.antigen_chain, dist_mat, indices, binding_interface)
        data[row.pdb] = res
    except Exception as e:
        print(f"Could not process {row.pdb}, {e}")

with open(data_path / "initial_data.pkl", 'wb') as fh:
    pickle.dump(data, fh)
