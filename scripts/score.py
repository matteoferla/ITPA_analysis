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
filename = 'hybrid_relaxed2.pdb'
params_filenames = ['ITP.params']
pose = pyrosetta.Pose()
if params_filenames:
    params_paths = pyrosetta.rosetta.utility.vector1_string()
    params_paths.extend(params_filenames)
    pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)

gnomad = '''p.Ala2Val
p.Leu5Ser
p.Lys8Arg
p.Lys9Arg
p.Val11Met
p.Ala17Ser
p.Leu20Met
p.Val24Ile
p.Pro32Thr
p.Thr34Ala
p.Thr34Asn
p.Leu35Ser
p.Val36Leu
p.Val36Glu
p.Lys39Asn
p.Pro43Leu
p.Glu48Gly
p.Pro49Leu
p.Glu51Lys
p.Ile54Val
p.Cys57Tyr
p.Val61Ile
p.Arg62Cys
p.Arg62Leu
p.Val64Leu
p.Gly66Arg
p.Val70Ile
p.Asn78Ser
p.Gly85Ala
p.Tyr87Asp
p.Tyr87Phe
p.Ile88Leu
p.Trp90Cys
p.Glu93Gly
p.Leu95Phe
p.His101Arg
p.Gln102His
p.Leu104Pro
p.Phe107Ser
p.Gly106Arg
p.Asp109Glu
p.Lys110Asn
p.Leu115Phe
p.Gly123Ser
p.Ser121Gly
p.Ser121Arg
p.Asp124Val
p.Pro125Thr
p.Pro125Ser
p.Gln127His
p.Pro128Leu
p.Val129Met
p.Val129Leu
p.Arg130Cys
p.Arg133Gly
p.Arg135Trp
p.Arg135Gln
p.Ser137Leu
p.Arg139Gln
p.Ile140Phe
p.Val141Met
p.Val141Gly
p.Ala142Thr
p.Ala142Val
p.Arg144Lys
p.Trp151Ter
p.Cys154Ala
p.Pro153His
p.Gln156His
p.Pro157Ala
p.Gly159Val
p.Tyr160Ser
p.Glu161Lys
p.Thr163Ala
p.Tyr164Ter
p.Ala165Thr
p.Ala170Thr
p.Ala170Val
p.Glu171Asp
p.Asn173Lys
p.Ala174Thr
p.Ala174Val
p.Val175Leu
p.Ser176Phe
p.Arg178Cys
p.Phe179Ile
p.Arg180Gln
p.Ala181Gly
p.Gln186Arg
p.Glu187Asp
p.Phe189Leu
p.Ser191Arg
p.Leu192Ser
p.Ala193Gly'''.split()

# ======================================================================================================================
from mutate import Variant

model = Variant(pose=pose) # .from_file(filename=pdb_filename)
model.scorefxn = pyrosetta.create_score_function('ref2015')
model.strict_about_starting_residue = True
model.score_mutations(sorted(gnomad, key=lambda v: int(v[5:-3])),
                        chain='A',
                        interfaces=(),
                        modelname='ITPA',
                        preminimise=False,
                        verbose=True)
print('DONE')