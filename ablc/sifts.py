#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import requests
import typing as t
from dataclasses import dataclass, field
from collections import defaultdict

ENTRY_ENDPOINT = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/"
SUMMARY_ENDPOINT = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/"


HEAVY_CHAIN_IDENTIFIERS = [
    "heavy",
    "vh",
    "nb",
    "body",
    "h-chain",
    "h chain",
    "chain h",
    "ig ",
    "ig-",
    "igg",
    "igm",
    "fab"
]

LIGHT_CHAIN_IDENTIFIERS = [
    "light",
    "l-chain",
    "l chain"
    "chain l",
    "kappa",
    "igk"
]


class SIFTSAntibodyError(Exception):
    def __init__(self, message):
        super().__init__(message)


def get_summary(pdb_id: str):
    response = requests.get(f"{SUMMARY_ENDPOINT}{pdb_id}")
    if response.status_code == 404:
        raise LookupError(f"No summary entry found in SIFTS for {pdb_id}")

    response = response.json()["pdb_id"]
    return response


@dataclass
class Entity:
    """
    SIFTS Entity
    """
    name: str
    tag: str
    chains: t.List[str]
    sequence: str = field(repr=False)
    pdb_sequence: str = field(repr=False)
    host: str
    host_tax_id: int = field(repr=False)
    org: str
    org_tax_id: int = field(repr=False)


def get_entities(pdb_id: str, nanobody: bool = False, require_antigen: bool = True) -> t.Dict[str, t.List[Entity]]:
    """
    Attempt to parse entities and their chain IDs from SIFTS

    Requires that all entities have a sequence record, ensuring protein antigens

    :param pdb_id: PDB ID
    :param nanobody: If true, does not attempt to find a light chain, default = False
    :param require_antigen: If true, does not attempt to find an antigen, default = False
    :return: A dictionary of ["heavy", "light", "antigen"] Entities
    :raises SIFTSAntibodyError: Could not find suitable SIFTS data
    """
    response = requests.get(f"{ENTRY_ENDPOINT}{pdb_id}")
    if response.status_code == 404:
        raise LookupError(f"No entity data found in SIFTS for {pdb_id}")

    entities = []
    found_hc = False
    found_lc = False
    found_ag = False
    response = response.json()[pdb_id]

    if nanobody:
        found_lc = True
    if not require_antigen:
        found_ag = True

    for entity in response:

        if "sequence" not in entity:
            continue
        try:
            ent = Entity(
                    name=entity["molecule_name"][0],
                    chains=entity["in_chains"],
                    tag="antigen",
                    sequence=entity["sequence"],
                    pdb_sequence=entity["pdb_sequence"],
                    host=entity["source"][0]["expression_host_scientific_name"],
                    host_tax_id=entity["source"][0]["expression_host_tax_id"],
                    org=entity["source"][0]["organism_scientific_name"],
                    org_tax_id=entity["source"][0]["tax_id"]
                )
        except (KeyError, IndexError):
            continue
        if nanobody:
            if any(s in ent.name.lower() for s in HEAVY_CHAIN_IDENTIFIERS):
                found_hc = True
                ent.tag = "heavy"
            if any(s in ent.org.lower() for s in ["camelid", "camelus", "lama", "vicugna"]):
                found_hc = True
                ent.tag = "heavy"
        else:
            if any(s in ent.name.lower() for s in HEAVY_CHAIN_IDENTIFIERS) \
                    and not any(s in ent.name.lower() for s in LIGHT_CHAIN_IDENTIFIERS):
                found_hc = True
                ent.tag = "heavy"

            elif any(s in ent.name.lower() for s in LIGHT_CHAIN_IDENTIFIERS):
                found_lc = True
                ent.tag = "light"

            elif "H" in ent.chains and not "L" in ent.chains:
                found_hc = True
                ent.tag = "heavy"
            elif "L" in ent.chains and not "H" in ent.chains:
                found_lc = True
                ent.tag = "light"

        entities.append(ent)

    if found_hc and found_lc and ((nanobody and len(entities) >= 2) or (not nanobody and len(entities) >= 3)):
        found_ag = True

    if found_ag and found_lc and found_hc:
        res = defaultdict(list)
        for entity in entities:
            res[entity.tag].append(entity)
        return res
    else:
        raise SIFTSAntibodyError(f"Could not find suitable SIFTS data for {pdb_id}, nanobody={nanobody}, "
                                 f"require_ag={require_antigen}, entities={entities}")
