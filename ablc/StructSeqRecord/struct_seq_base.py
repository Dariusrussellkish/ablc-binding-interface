#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>
import os
from abc import ABC
from enum import Enum
import typing as t
from ..sifts import Entity as SIFTSEntity
from ..anarci import ANARCI

from Bio.SeqIO import SeqRecord
from Bio.PDB.Residue import Residue


class CHAINS(Enum):
    LIGHT = "light"
    HEAVY = "heavy"
    ANTIGEN = "antigen"


class BaseStructSeqRecord(ABC):
    def __init__(self,
                 pdb_id: str,):
        self.pdb_id: str = pdb_id


class BaseUnsequencedEntityRecord(ABC):
    def __init__(self,
                 chain_names: t.Set[str],
                 sifts_records: t.List[SIFTSEntity]
                 ):
        self.chain_names = chain_names
        self.sifts_records = sifts_records


class BaseSequencedEntityRecord(BaseUnsequencedEntityRecord):
    def __init__(self,
                 chain_names: t.Set[str],
                 sifts_records: t.List[SIFTSEntity],
                 sequence: SeqRecord,
                 residues: t.Sequence[Residue],
                 anarci: ANARCI
                 ):
        super().__init__(chain_names, sifts_records)
        self.sequence = sequence
        self.anarci = anarci
        self.residues = residues
