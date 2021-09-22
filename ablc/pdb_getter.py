#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>
import os
from urllib.request import urlopen
from urllib.error import HTTPError
from shutil import copyfileobj
from tempfile import NamedTemporaryFile
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.Structure import Structure
from pathlib import Path
from os import path, makedirs
import gzip

pdb_parser = PDBParser(QUIET=True)

BIOASSEBLY_BASE_URL = "https://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/all"


class PDBGetterError(Exception):
    def __init__(self, msg):
        super().__init__(msg)


def ensure_dir(pth: os.PathLike[str]) -> None:
    makedirs(Path(pth).resolve().absolute(), exist_ok=True)


def bioassembly_template(pdb: str, assembly_id=1) -> str:
    return f"{BIOASSEBLY_BASE_URL}/{pdb}.pdb{assembly_id}.gz"


def _get_bioassembly_pdb(pdb: str,
                        assembly_id = 1,
                        store_location: os.PathLike = None,
                        retry = False) -> Structure:
    try:
        if store_location is None:
            with urlopen(bioassembly_template(pdb, assembly_id=assembly_id)) as fsrc, NamedTemporaryFile() as fdst:
                with gzip.open(fsrc, 'rb') as fin:
                    copyfileobj(fin, fdst)

                structure = pdb_parser.get_structure(pdb, fdst.name)
                return structure
        else:
            ensure_dir(store_location)
            with urlopen(bioassembly_template(pdb, assembly_id=assembly_id)) as fsrc, \
                    open(path.join(store_location, f"{pdb}_{assembly_id}.pdb"), "wb") as fdst:
                with gzip.open(fsrc, 'rb') as fin:
                    copyfileobj(fin, fdst)

                structure = pdb_parser.get_structure(pdb, fdst.name)
                return structure
    except HTTPError as e:
        if retry:
            raise PDBGetterError(f"No parseable bioassembly for {pdb} found")
        else:
            raise PDBGetterError(f"Bioassembly {assembly_id} for {pdb} not found")
    except PDBConstructionException as e:
        _get_bioassembly_pdb(pdb, assembly_id=assembly_id + 1, store_location=store_location, retry=True)


def get_bioassembly_pdb(pdb: str,
                        assembly_id = 1,
                        store_location: os.PathLike = None) -> Structure:
    """
    Attempt to retrieve a bioassembly PDB file from the PDB

    :param pdb: PDB ID
    :param assembly_id: Optional, assembly ID, default = 1
    :param store_location: Optional, store physical PDB files, default = None
    :raises PDBGetterError: No (valid) bioassembly exists for PDB ID
    :return: Biopython PDB Structure
    """
    return _get_bioassembly_pdb(pdb,
                                assembly_id=assembly_id,
                                store_location=store_location)
