# ED without membrane
# ======== init pyrosetta ==============================================================================================
import pyrosetta
from init import make_option_string


pyrosetta.init(extra_options=make_option_string(no_optH=False,
                                                ex1=None,
                                                ex2=None,
                                                mute='all',
                                                ignore_unrecognized_res=True,
                                                load_PDB_components=False,
                                                ignore_waters=False)
               )

# ======== Pose ========================================================================================================
filename = ' hybrid.map_aligned.pdb'
params_filenames = ['ITP.params']
pose = pyrosetta.Pose()
if params_filenames:
    params_paths = pyrosetta.rosetta.utility.vector1_string()
    params_paths.extend(params_filenames)
    pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)

# ======== Check ========================================================================================================
import nglview

view = nglview.show_rosetta(pose)
view.add_representation('spacefill', selection='ion')
view.add_representation('hyperball', selection='water')
view.add_representation('hyperball', selection='ITT', color='white')
view

# ======== scorefxn ====================================================================================================

scorefxnED = pyrosetta.get_fa_scorefxn()  # ref2015 or franklin2019?


def get_ED(pose: pyrosetta.Pose, map_filename: str):
    ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_filename)  # from PDBe
    sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
    sdsm.apply(pose)
    return ED


def add_ED(ED, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction):
    elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
    scorefxn.set_weight(elec_dens_fast, 30)
    return scorefxn


# ref2015 or franklin2019?
ED = get_ED(pose, '2j4e.ccp4')
scorefxnED = add_ED(ED, pyrosetta.create_score_function('ref2015'))
scorefxnED_cart = add_ED(ED, pyrosetta.create_score_function('ref2015_cart'))

# ======== relax =======================================================================================================
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED_cart, 5)
relax.minimize_bond_angles(True)
relax.minimize_bond_lengths(True)
relax.apply(pose)

for w in (30, 20, 10):
    elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
    scorefxnED.set_weight(elec_dens_fast, w)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 5)
    #     relax.minimize_bond_angles(True)
    #     relax.minimize_bond_lengths(True)
    relax.apply(pose)
    print('Match to ED', ED.matchPose(pose))
    print('ED weight', scorefxnED.weights()[elec_dens_fast])
    print('score', scorefxnED(pose))
reg_scorefxn = pyrosetta.get_fa_scorefxn()
print('How good is the density match:', ED.matchPose(pose))
print(reg_scorefxn.get_name(), reg_scorefxn(pose))

# ======== final =======================================================================================================
scorefxn = pyrosetta.create_score_function('ref2015')
print(scorefxn(pose))
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
relax.apply(pose)
print(scorefxn(pose))
pose.dump_pdb('hybrid_relaxed.pdb')