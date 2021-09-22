#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import pathlib
import pickle
import os
import traceback


from ablc.StructSeqRecord import *
from ablc.structure_analysis import *
from ablc.df_query import *
from ablc.sifts import SIFTSAntibodyError

from Bio.PDB.StructureBuilder import PDBConstructionWarning
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO
from Bio import AlignIO
import warnings
from joblib import Parallel, delayed

warnings.simplefilter('ignore', category=PDBConstructionWarning)

current_dir = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
data_path = current_dir / "data"

df = pd.read_csv(data_path / "sabdab_nr_h_l_bound_summary.tsv", sep='\t')

df = filter_for_unique_chains(df)
df = filter_for_peptide_antigens(df)
df = combine_dup_rows(df)

data = dict()
rows = df.shape[0]

good_records = {}
obsolete_records = {}
mismatched_records = {}
failed_seq_records = {}
all_records = [good_records, obsolete_records, mismatched_records, failed_seq_records]


def process_pdb(row):
    try:
        record = build_struct_seq_record(row.pdb,
                                         model=row.model,
                                         default_save_dir=(data_path / "pdb"))
        return row.pdb, record
    except Exception as e:
        print(f"Could not process {row.pdb}, {e}\n")
        if not isinstance(e, SIFTSAntibodyError):
            traceback.print_exc()
            print("\n\n")
        return row.pdb, None


if not os.path.exists('parsed_records.pkl'):
    res = Parallel(n_jobs=8)(delayed(process_pdb)(row) for row in df.itertuples(index=False))
    for pdb, record in res:
        if record is None:
            failed_seq_records[pdb] = record
        if isinstance(record, CompleteStructSeqRecord):
            good_records[pdb] = record
        elif isinstance(record, MismatchStructSeqRecord):
            mismatched_records[pdb] = record
        elif isinstance(record, MissingSeqStructSeqRecord):
            failed_seq_records[pdb] = record
        elif isinstance(record, ObsoleteStructSeqRecord):
            obsolete_records[pdb] = record

    print(f"Good: {len(good_records)} | Mismatch: {len(mismatched_records)} | Failed: {len(failed_seq_records)} "
          f"| Obsolete: {len(obsolete_records)} | Total: {sum([len(x) for x in all_records])}")

    with open('parsed_records.pkl', 'wb') as fh:
        pickle.dump(good_records, fh)
else:
    with open('parsed_records.pkl', 'rb') as fh:
        good_records = pickle.load(fh)


def process_record(record: CompleteStructSeqRecord):
    light_chains = {record.pdb_struct[record.model][c] for c in record.entity_records["light"].chain_names}
    antigen_chains = {record.pdb_struct[record.model][c] for c in record.entity_records["antigen"].chain_names}
    record.binding_interfaces["light-antigen"] = get_binding_partners(light_chains, antigen_chains, cutoff=6.0)
    return record


if not os.path.exists('processed_records.pkl'):
    processed_records = Parallel(n_jobs=8)(delayed(process_record)(record) for record in good_records.values())
    for record in processed_records:
        good_records[record.pdb_id] = record
    with open('processed_records.pkl', 'wb') as fh:
        pickle.dump(good_records, fh)
else:
    with open('processed_records.pkl', 'rb') as fh:
        good_records = pickle.load(fh)

with open('unaligned_good_records.fasta', 'w') as fh:
    SeqIO.write([record.entity_records["light"].sequence for record in good_records.values()], fh, format='fasta')

clustalo_align = ClustalOmegaCommandline(infile="unaligned_good_records.fasta",
                                         outfile="aligned_good_records.fasta",
                                         force=True,
                                         verbose=True, auto=True)
clustalo_align()
alignment = AlignIO.read("aligned_good_records.fasta", "fasta")
