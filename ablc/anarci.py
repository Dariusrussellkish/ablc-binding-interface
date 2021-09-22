#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

from enum import Enum
from dataclasses import dataclass
import typing as t
import subprocess


class IMGT(Enum):
    CDR1START = "27"
    CDR1END = "38"
    CDR2START = "56"
    CDR2END = "65"
    CDR3START = "105"
    CDR3END = "117"


@dataclass
class ANARCI:
    species: str
    chain: str
    e_val: float
    score: float
    start: int
    end: int
    alignment: t.Dict[str, str]
    anarci_to_sequence_idx: t.Dict[str, int]
    CDR1: t.Dict[int, str]
    CDR2: t.Dict[int, str]
    CDR3: t.Dict[int, str]
    CDR1_to_alignment: t.Dict[int, str]
    CDR2_to_alignment: t.Dict[int, str]
    CDR3_to_alignment: t.Dict[int, str]


class ANARCIError(Exception):
    def __init__(self, message):
        super().__init__(message)


def single_anarci_alignment(seq: str):
    res = subprocess.run(
        f"ANARCI -i {seq} -s imgt",
        shell=True,
        capture_output=True,
        text=True,
    )

    res.check_returncode()

    species = None
    chain = None
    e_val = None
    score = None
    start = None
    end = None
    alignment = dict()
    anarci_to_sequence_idx = dict()
    cdr1 = dict()
    cdr2 = dict()
    cdr3 = dict()
    cdr1_to_alignment = dict()
    cdr2_to_alignment = dict()
    cdr3_to_alignment = dict()

    for i, line in enumerate(res.stdout.splitlines()):
        if i < 5 or i == 6:
            continue
        line = line.strip()
        if line == "//":
            break
        if "#|" in line:
            _, species, chain, e_val, score, start, end, _ = line.split("|")
            e_val = float(e_val)
            score = float(score)
            start = int(start)
            end = int(end)
        else:
            vals = line.split()
            if len(vals) == 3:
                _, imgt, res = vals
                alignment[imgt] = res
            elif len(vals) == 4:
                _, imgt, sub, res = vals
                alignment[imgt+sub] = res
            else:
                raise ANARCIError(f"Unexpected number of values in line: {i+1} of {res.stdout} for sequence: {seq}")

    if score is None:
        raise ANARCIError(f"ANARCI unsuccessful: {res.stdout}")

    phase = -1
    gap_offset = 0
    for i, (imgt, res) in enumerate(alignment.items()):
        if imgt == IMGT.CDR1START.value:
            phase = 1
        elif imgt == IMGT.CDR1END.value:
            phase = -1
        elif imgt == IMGT.CDR2START.value:
            phase = 2
        elif imgt == IMGT.CDR2END.value:
            phase = -1
        elif imgt == IMGT.CDR3START.value:
            phase = 3
        elif imgt == IMGT.CDR3END.value:
            phase = -1

        if res == "-":
            gap_offset += 1
            continue

        anarci_to_sequence_idx[imgt] = start + i - gap_offset

        if phase == 1:
            cdr1[start + i - gap_offset] = res
            cdr1_to_alignment[start + i - gap_offset] = imgt
        elif phase == 2:
            cdr2[start + i - gap_offset] = res
            cdr2_to_alignment[start + i - gap_offset] = imgt
        elif phase == 3:
            cdr3[start + i - gap_offset] = res
            cdr3_to_alignment[start + i - gap_offset] = imgt

    return ANARCI(
        species,
        chain,
        e_val,
        score,
        start,
        end,
        alignment,
        anarci_to_sequence_idx,
        cdr1,
        cdr2,
        cdr3,
        cdr1_to_alignment,
        cdr2_to_alignment,
        cdr3_to_alignment
    )