#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import typing as t
from contextlib import redirect_stdout
from io import StringIO
from os import PathLike

import requests
from Bio.PDB import PDBList, PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import is_aa, three_to_one
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


from .struct_seq_base import BaseStructSeqRecord, BaseUnsequencedEntityRecord
from .struct_seq_complete import CompleteStructSeqRecord, AntigenChainEntityRecord, LightChainEntityRecord, \
    HeavyChainEntityRecord
from .struct_seq_partials import MismatchStructSeqRecord, ObsoleteStructSeqRecord, MissingSeqStructSeqRecord
from .struct_seq_partials import NonExistentPDBException
from .binding_interface import BindingInterface
from ..sifts import get_entities
from ..anarci import single_anarci_alignment
from ..structure_analysis import get_binding_partners
from ..pdb_getter import get_bioassembly_pdb, PDBGetterError

pdb_parser = PDBParser(QUIET=True)
pdb_list = PDBList(verbose=False)
obsolete_pdbs = [s.lower() for s in requests.get("https://data.rcsb.org/rest/v1/holdings/removed/entry_ids").json()]


def download_or_open_pdb(pdb_id: str,
                         default_save_dir: t.Union[str, PathLike] = "pdb") -> t.Optional[Structure]:
    """
    Fetch PDB (physically download it) and return a Biopython PDB struct

    :param pdb_id: PDB ID
    :param default_save_dir: directory to save physical PDB files
    :return: Optional biopython structure
    """
    with redirect_stdout(StringIO()) as _:
        if pdb_id in obsolete_pdbs:
            return None
        pdb_location = pdb_list.retrieve_pdb_file(pdb_id, pdir=default_save_dir, file_format="pdb")
        structure = pdb_parser.get_structure(pdb_id, pdb_location)
        if structure is None:
            raise NonExistentPDBException(f"Could not find PDB file for id: {pdb_id}")
        return structure


def get_residues(chain_names: t.Set[str], model: Model) -> t.Iterable[Residue]:
    for chain in chain_names:
        for residue in model[chain]:
            if is_aa(residue.get_resname(), standard=True):
                yield residue
            else:
                continue


class MismatchAlignmentError(Exception):
    def __init__(self, message):
        super().__init__(message)


def build_struct_seq_record(pdb_id: str,
                            model=0,
                            default_save_dir: t.Union[str, PathLike] = "pdb",
                            nanobody=False,
                            require_antigen=True,
                            ) -> t.Union[CompleteStructSeqRecord,
                                         MismatchStructSeqRecord,
                                         ObsoleteStructSeqRecord,
                                         MissingSeqStructSeqRecord]:
    """
    Create a Structure-Sequence Record for a PDB ID

    :param pdb_id: PDB ID
    :param model: Optional, model index of the PDB struct, default = 0
    :param default_save_dir: Optional, directory to save physical PDB files, default = "pdb"
    :param nanobody: Optional, if true, don't look for light chains, default = False
    :return: Structure-Sequence Record
    """

    if pdb_id in obsolete_pdbs:
        return ObsoleteStructSeqRecord(pdb_id)

    # try to get a SIFTS record for the PDB ID
    sifts_records = get_entities(pdb_id, nanobody=nanobody, require_antigen=require_antigen)

    # Here we are trying to get a valid bioassembly PDB that matches the SIFTS record.
    # sometimes they're malformed so we might as well repeat for all available bioassemblies
    successful_matching = False
    assembly_id = 1
    while not successful_matching:
        try:
            struct = get_bioassembly_pdb(pdb_id, assembly_id, store_location=default_save_dir)
            smodel = struct[model]
        except PDBGetterError as e:
            raise e

        # Mapping SIFTS chain records to the physical PDB structure
        pdb_ags = set()
        sifts_ags = set([c for r in sifts_records["antigen"] for c in r.chains])
        pdb_hcs = set()
        sifts_hcs = set([c for r in sifts_records["heavy"] for c in r.chains])
        pdb_lcs = set()
        if not nanobody:
            sifts_lcs = set([c for r in sifts_records["light"] for c in r.chains])
        else:
            sifts_lcs = set()
        for chain in set([c.get_id() for c in smodel.get_chains()]):
            if chain in sifts_ags:
                pdb_ags.add(chain)
            elif chain in sifts_hcs:
                pdb_hcs.add(chain)
            elif not nanobody and chain in sifts_lcs:
                pdb_lcs.add(chain)
        try:
            # in the case of dimer/trimer bioassemblies we only care about one of the HCs
            pdb_hcs = set(pdb_hcs.pop())
            successful_matching = True
        except KeyError:
            assembly_id += 1

    entity_dict = dict()
    if require_antigen:
        entity_dict["antigen"] = AntigenChainEntityRecord(pdb_ags, sifts_records["antigen"])

    # since we took only one HC, we need to find its matching pair
    if not nanobody and len(pdb_lcs) != 1:
        lc_hc_contacts = {}
        for lc in pdb_lcs:
            lc_hc_contacts[lc] = len(get_binding_partners({smodel[lc]}, {smodel[c] for c in pdb_hcs}).interface_pairs)
        lc = max(lc_hc_contacts, key=lc_hc_contacts.get)
        if lc_hc_contacts[lc] < 1:
            raise MismatchAlignmentError(f"Could not find an appropriate lc, hc pair for pdb {pdb_id}")

        pdb_lcs = set(lc)

    hc_seq = "".join(three_to_one(r.get_resname()) for r in get_residues(pdb_hcs, smodel))
    hc_anarci = single_anarci_alignment(hc_seq)
    entity_dict["heavy"] = HeavyChainEntityRecord(pdb_hcs,
                                                  sifts_records["heavy"],
                                                  SeqRecord(seq=Seq(hc_seq),
                                                            id=f"{pdb_id}|Heavy|{','.join(pdb_hcs)}|{len(hc_seq)}",
                                                            description=""),
                                                  list(get_residues(pdb_hcs, smodel)),
                                                  hc_anarci)
    if not nanobody:
        lc_seq = "".join("".join(three_to_one(r.get_resname()) for r in get_residues(pdb_lcs, smodel)))
        lc_anarci = single_anarci_alignment(lc_seq)
        entity_dict["light"] = LightChainEntityRecord(pdb_lcs,
                                                      sifts_records["light"],
                                                      SeqRecord(seq=Seq(lc_seq),
                                                                id=f"{pdb_id}|Light|{','.join(pdb_lcs)}|{len(lc_seq)}",
                                                                description=""),
                                                      list(get_residues(pdb_lcs, smodel)),
                                                      lc_anarci)

    if require_antigen and len(pdb_ags) == 0:
        return MissingSeqStructSeqRecord(pdb_id,
                                         struct,
                                         "antigen")

    return CompleteStructSeqRecord(
        pdb_id,
        struct,
        model,
        sifts_records,
        entity_dict,
    )
