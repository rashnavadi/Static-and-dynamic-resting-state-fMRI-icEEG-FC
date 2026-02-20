#!/usr/bin/env python3
import os
import re
import glob
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import plotting


# =========================
# USER SETTINGS
# =========================
subj = "ICE013"

base_dir = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}"
preproc_dir = os.path.join(base_dir, "4_MRI", "2_Functionals", "2_PreProcessing")

anat_native_nii = os.path.join(
    base_dir, "4_MRI", "1_Anatomicals", "3D_Anat", f"{subj}_3Danat_brain.nii.gz"
)

coords_xlsx = f"/Volumes/MIND/ICE/original_ICE/MNI_coordinates/{subj}_channel_info.xlsx"

out_dir = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_outputs_py"
os.makedirs(out_dir, exist_ok=True)


# =========================
# HELPERS
# =========================
def find_reg_dir(preproc_dir):
    # Find Run*.feat, exclude anything containing "Exclude"
    run_feats = sorted(glob.glob(os.path.join(preproc_dir, "Run*.feat")))
    run_feats = [r for r in run_feats if "Exclude" not in os.path.basename(r)]
    if not run_feats:
        raise FileNotFoundError(f"No non-Exclude Run*.feat found in: {preproc_dir}")
    # pick the first usable run
    reg_dir = os.path.join(run_feats[0], "reg")
    if not os.path.isdir(reg_dir):
        raise FileNotFoundError(f"reg folder not found: {reg_dir}")
    return reg_dir, os.path.basename(run_feats[0])


def guess_xyz_cols(df):
    cols = list(df.columns)
    def n(s): return re.sub(r"[^a-z0-9]+", "", str(s).lower())

    nmap = {c: n(c) for c in cols}

    # Prefer explicit MNI columns if they exist
    def pick(axis):
        prefs = [f"mni{axis}", f"mni_{axis}", axis]
        prefs = [n(p) for p in prefs]
        for c in cols:
            if nmap[c] in prefs:
                return c
        # fallback: any column containing axis and looking coordinate-like
        for c in cols:
            s = nmap[c]
            if axis in s and ("mni" in s or "coord" in s or "mm" in s):
                return c
        return None

    cx, cy, cz = pick("x"), pick("y"), pick("z")
    if cx is None or cy is None or cz is None:
        raise ValueError(f"Could not detect x/y/z columns. Columns: {cols}")
    return cx, cy, cz


def run_std2imgcoord_to_mm(mni_xyz_mm, reg_dir):
    """
    Convert MNI(mm) coords -> highres(mm) coords using:
      std2imgcoord -img highres -std standard -xfm highres2standard.mat -mm
    """
    highres = os.path.join(reg_dir, "highres.nii.gz")
    standard = os.path.join(reg_dir, "standard.nii.gz")
    xfm = os.path.join(reg_dir, "highres2standard.mat")  # IMPORTANT direction

    for p in [highres, standard, xfm]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required reg file: {p}")

    tmp_in = os.path.join(out_dir, "tmp_mni_coords.txt")
    tmp_out = os.path.join(out_dir, "tmp_highres_mm.txt")
    np.savetxt(tmp_in, mni_xyz_mm, fmt="%.6f")

    cmd = [
        "std2imgcoord",
        "-img", highres,
        "-std", standard,
        "-xfm", xfm,
        "-mm"
    ]

    # std2imgcoord reads coords from stdin
    with open(tmp_in, "r") as fin, open(tmp_out, "w") as fout:
        p = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=subprocess.PIPE, text=True)

    if p.returncode != 0:
        raise RuntimeError(f"std2imgcoord failed:\n{p.stderr}")

    xyz_hi_mm = np.loadtxt(tmp_out)
    if xyz_hi_mm.ndim == 1:
        xyz_hi_mm = xyz_hi_mm.reshape(1, 3)
    if xyz_hi_mm.shape[1] != 3:
        raise ValueError(f"Unexpected output shape from std2imgcoord: {xyz_hi_mm.shape}")

    return xyz_hi_mm, highres


def affines_close(A, B, tol_mm=2.0):
    # rough check: max absolute difference
    return np.max(np.abs(A - B)) < tol_mm


# =========================
# MAIN
# =========================
reg_dir, run_used = find_reg_dir(preproc_dir)
print(f"[INFO] Using reg dir: {reg_dir} (from {run_used})")

df = pd.read_excel(coords_xlsx)
cx, cy, cz = guess_xyz_cols(df)

mni_xyz = df[[cx, cy, cz]].to_numpy(dtype=float)
good = np.all(np.isfinite(mni_xyz), axis=1)
df_use = df.loc[good].copy()
mni_xyz = mni_xyz[good]
print(f"[INFO] Loaded {len(mni_xyz)} finite MNI coords from Excel.")

# Convert MNI(mm) -> highres(mm) using FSL
hi_xyz_mm, highres_nii = run_std2imgcoord_to_mm(mni_xyz, reg_dir)
print("[OK] Converted MNI -> highres (mm) using std2imgcoord + highres2standard.mat")

# Save CSV
out_csv = os.path.join(out_dir, f"{subj}_MNI_to_highres_mm.csv")
out = df_use.copy()
out["mni_x_mm"] = mni_xyz[:, 0]
out["mni_y_mm"] = mni_xyz[:, 1]
out["mni_z_mm"] = mni_xyz[:, 2]
out["highres_x_mm"] = hi_xyz_mm[:, 0]
out["highres_y_mm"] = hi_xyz_mm[:, 1]
out["highres_z_mm"] = hi_xyz_mm[:, 2]
out.to_csv(out_csv, index=False)
print("[OK] Saved:", out_csv)

# QC 1: overlay on highres (this is the most trustworthy QC)
qc1 = os.path.join(out_dir, f"{subj}_QC_on_highres.png")
disp = plotting.plot_anat(highres_nii, title=f"{subj}: electrodes on FEAT highres (MNI→highres via std2imgcoord)")
disp.add_markers(hi_xyz_mm, marker_size=40)
disp.savefig(qc1, dpi=200)
disp.close()
print("[OK] Saved:", qc1)

# QC 2: overlay on your native 3Danat (only if aligned)
anat_img = nib.load(anat_native_nii)
hi_img = nib.load(highres_nii)

qc2 = os.path.join(out_dir, f"{subj}_QC_on_native3Danat.png")
if affines_close(anat_img.affine, hi_img.affine, tol_mm=5.0):
    disp = plotting.plot_anat(anat_native_nii, title=f"{subj}: electrodes on native 3Danat (highres assumed aligned)")
    disp.add_markers(hi_xyz_mm, marker_size=40)
    disp.savefig(qc2, dpi=200)
    disp.close()
    print("[OK] Saved:", qc2)
else:
    print("[WARN] highres.nii.gz affine != native 3Danat affine (not trivially aligned).")
    print("       Use the highres QC as ground truth, or we need a highres→3Danat mapping step.")
    print("       (I can add that next.)")

print("[DONE]")
