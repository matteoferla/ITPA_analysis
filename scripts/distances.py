## code from a Jupyter notebook

import re
import pyrosetta
from pyrosetta_help import make_option_string
pyrosetta.init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
               )

# ======== Pose ========================================================================================================
filename = 'model.pdb'
params_filenames = ['ITP.params']
pose = pyrosetta.Pose()
if params_filenames:
    params_paths = pyrosetta.rosetta.utility.vector1_string()
    params_paths.extend(params_filenames)
    pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)

pathogenic = ['👾👾👾', '👾👾👾', '👾👾👾']


pdb2pose = pose.pdb_info().pdb2pose


def distance(atom_A: pyrosetta.rosetta.core.conformation.Atom,
             atom_B: pyrosetta.rosetta.core.conformation.Atom) -> float:
    assert isinstance(atom_A, pyrosetta.rosetta.core.conformation.Atom)
    assert isinstance(atom_B, pyrosetta.rosetta.core.conformation.Atom)
    bond = atom_A.xyz() - atom_B.xyz()
    return bond.norm()


itt_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector('ITT')
mg_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector('MG')
or_sele = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(itt_sele, mg_sele)
ligs = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(or_sele.apply(pose))

for mutant in pathogenic:
    print(mutant, flush=True)
    r = pdb2pose(chain='A', res=int(re.search('(\d+)', mutant).group(1)))
    assert r, f'mutant {mutant} ??'
    residue = pose.residue(r)
    for lig_i in ligs:
        lig = pose.residue(lig_i)
        distances = []
        for a in range(1, residue.natoms() + 1):
            for b in range(1, lig.natoms() + 1):
                distances.append(distance(residue.atom(a), lig.atom(b)))
        print(mutant, lig.name3(), min(distances))
