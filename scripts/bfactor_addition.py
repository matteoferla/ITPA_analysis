import pyrosetta
from typing import Dict

def add_bfactor_from_consurf(pose: pyrosetta.Pose, grades_filename: str) -> Dict[str, float]:
    """
    Adds the bfactors from a consurf run to a pose based on the PDB residues.
    The file ``consurf.grades`` is parsed weirdly as it is a weird file.
    ``replace_res_remap_bfactors`` or ``set_bfactors_from_atom_id_map``
    were not used but may have been a cleaner strategy. This was quicker to write.
    """
    # parse file
    consurf_scores = {}
    with open(grades_filename) as r:
        for row in r:
            row = row.strip()
            if not len(row) or not row[0].isdigit():
                continue
            parts = row.split()
            consurf_scores[parts[2]] = float(parts[3])
    # add to pose
    pdb_info = pose.pdb_info()
    pdb2pose = pdb_info.pdb2pose
    for con_res in consurf_scores:
        pose_res = pdb2pose(res=int(con_res[3:-2]), chain=con_res[-1])
        assert pose.residue(pose_res).name3() == con_res[:3], f'{pose.residue(pose_res).name3()} ≠ {con_res[:3]}'
        for i in range(pose.residue(pose_res).natoms()):
            pdb_info.bfactor(pose_res, i, consurf_scores[con_res])
    return consurf_scores
