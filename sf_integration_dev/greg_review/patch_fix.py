# Create a fixed copy: add a POST-GROWTH density-cap re-application so end-of-step
# SDI cannot exceed SDIMAX (iterate-to-convergence in one extra pass). Works on the
# COPY in track_site; originals untouched.
import re
src="/fs/scratch/PUOM0008/crsfaaron/track_site/constrained_projection.py"
dst="/fs/scratch/PUOM0008/crsfaaron/track_site/constrained_projection_fixed.py"
code=open(src).read()
# Inject right before 'traj.append(snapshot(yr))'. Re-cap TPA after growth+survival.
inject = (
'        # ---- (4b) FIX (Track C): POST-GROWTH density-cap re-application ----\n'
'        # After DBH growth + survival, QMD has risen so stand SDI can drift back\n'
'        # above SDIMAX. Re-apply the Reineke cap on the grown state so END-OF-STEP\n'
'        # SDI <= SDIMAX (removes the ~4% within-step rescale slack). One pass is\n'
'        # enough because re-scaling TPA does not change QMD (SDI is linear in TPA).\n'
'        if constrained:\n'
'            sdi_post = _imperial_sdi(dbh, tpa)          # (D,) grown-state imperial SDI\n'
'            for d in range(n_draws):\n'
'                s_d = float(sdi_post[d])\n'
'                if s_d > sdimax[d] and s_d > 0:\n'
'                    tpa[d] *= sdimax[d] / s_d\n'
'\n'
'        traj.append(snapshot(yr))'
)
assert code.count('        traj.append(snapshot(yr))')==1, "anchor count!="+str(code.count('        traj.append(snapshot(yr))'))
code=code.replace('        traj.append(snapshot(yr))', inject)
open(dst,"w").write(code)
print("wrote",dst)
