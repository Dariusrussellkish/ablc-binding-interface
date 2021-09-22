#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

from .struct_seq_base import BaseStructSeqRecord
from Bio.PDB import Structure


class ObsoleteStructSeqRecord(BaseStructSeqRecord):
    def __init__(self,
                 pdb_name: str,
                 ):
        super().__init__(pdb_name)


class MissingSeqStructSeqRecord(BaseStructSeqRecord):
    def __init__(self,
                 pdb_name: str,
                 pdb_struct: Structure,
                 failed_sequence: str
                 ):
        super().__init__(pdb_name)
        self.failed_sequence = failed_sequence
        self.pdb_struct = pdb_struct

    def __repr__(self):
        return f"MissingSeqStructSeqRecord<pdb: {self.pdb_id}, failed: {self.failed_sequence}>"


class MismatchStructSeqRecord(BaseStructSeqRecord):
    def __init__(self,
                 pdb_name: str,
                 pdb_struct: Structure,
                 failed_struct: str,
                 ):
        super().__init__(pdb_name)
        self.pdb_struct = pdb_struct
        self.failed_struct = failed_struct


class NonExistentPDBException(Exception):
    def __init__(self, message):
        super().__init__(message)
