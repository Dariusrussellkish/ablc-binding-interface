#  Copyright (c) 2021, Darius Russell Kish <dariusrussellkish@g.harvard.edu>

import xmlrpc.client as xmlrpclib


from .StructSeqRecord import CompleteStructSeqRecord


def show_bound_residues(record: CompleteStructSeqRecord):
    cmd = xmlrpclib.ServerProxy('http://localhost:9123')

    cmd.delete("all")
    cmd.fetch(record.pdb_id)
    cmd.hide("all")
    lc_set = set([p[0] for p in record.binding_interfaces["heavy-light"].interface_pairs])
    ag_set = set([p[1] for p in record.binding_interfaces["heavy-light"].interface_pairs])

    cmd.show("cartoon", " or ".join([f"chain {chain}" for chain in record.entity_records["light"].chain_names]))
    cmd.color("green", f"name CA and ({' or '.join([f'chain {chain}' for chain in record.entity_records['light'].chain_names])})")
    cmd.show("cartoon", " or ".join([f"chain {chain}" for chain in record.entity_records["antigen"].chain_names]))
    cmd.color("red", f"name CA and ({' or '.join([f'chain {chain}' for chain in record.entity_records['antigen'].chain_names])})")

    resid_selections = []

    for resid in lc_set:
        sele = f"(chain {resid.get_parent().id} and resid {resid.id[1]})"
        resid_selections.append(sele)
        cmd.show("sticks", sele)
        cmd.label2(f"name CB and {sele}", "resn")
        cmd.color("cyan", f"elem c and {sele}")
        cmd.color("red", f"elem o and {sele}")
        cmd.color("blue", f"elem n and {sele}")

    for resid in ag_set:
        sele = f"(chain {resid.get_parent().id} and resid {resid.id[1]})"
        resid_selections.append(sele)
        cmd.show("sticks", sele)
        cmd.label2(f"name CB and {sele}", "resn")
        cmd.color("pink", f"elem c and {sele}")
        cmd.color("red", f"elem o and {sele}")
        cmd.color("blue", f"elem n and {sele}")

    cmd.color("green", f"name CA and ({' or '.join([f'chain {chain}' for chain in record.entity_records['light'].chain_names])})")
    cmd.color("red", f"name CA and ({' or '.join([f'chain {chain}' for chain in record.entity_records['antigen'].chain_names])})")

    cmd.orient(" or ".join([f"{sele}" for sele in resid_selections]))

    _ = input("Done?> ")
    if _.strip():
        return "Done"
    return None
