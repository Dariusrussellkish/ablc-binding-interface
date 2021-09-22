#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

from .struct_seq_base import BaseStructSeqRecord, BaseUnsequencedEntityRecord, BaseSequencedEntityRecord
from ..sifts import Entity as SIFTSEntity
from ..anarci import ANARCI
from .binding_interface import BindingInterface

import typing as t

from Bio.PDB import Structure
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Residue import Residue


class LightChainEntityRecord(BaseSequencedEntityRecord):
    def __init__(self,
                 chain_names: set[str],
                 sifts_records: list[SIFTSEntity],
                 sequence: SeqRecord,
                 residues: t.Sequence[Residue],
                 anarci: ANARCI
                 ):
        super().__init__(chain_names, sifts_records, sequence, residues, anarci)


class HeavyChainEntityRecord(BaseSequencedEntityRecord):
    def __init__(self,
                 chain_names: set[str],
                 sifts_records: list[SIFTSEntity],
                 sequence: SeqRecord,
                 residues: t.Sequence[Residue],
                 anarci: ANARCI
                 ):
        super().__init__(chain_names, sifts_records, sequence, residues, anarci)


class AntigenChainEntityRecord(BaseUnsequencedEntityRecord):
    def __init__(self,
                 chain_names: set[str],
                 sifts_records: list[SIFTSEntity],
                 ):
        super().__init__(chain_names, sifts_records)


class CompleteStructSeqRecord(BaseStructSeqRecord):
    def __init__(self,
                 pdb_name: str,
                 pdb_struct: Structure,
                 model: int,
                 sifts_records: t.Dict[str, t.List[SIFTSEntity]],
                 entity_records: t.Dict[str, t.Union[LightChainEntityRecord, HeavyChainEntityRecord, AntigenChainEntityRecord]]
                 ):
        super().__init__(pdb_name)
        self.pdb_struct = pdb_struct
        self.model = model
        self.sifts_entity = sifts_records
        self.binding_interfaces: t.Dict[str, BindingInterface] = {}
        self.entity_records = entity_records

