#!/usr/bin/env python3
import os
import re
import glob
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import plotting
import matplotlib.pyplot as plt


# convert Excel MNI coords → native by using FEAT’s reg/ folder and std2imgcoord to map MNI(mm) 
# coordinates into the subject’s highres anatomical mm space, then QC by plotting the points on highres.nii.gz.
# =========================
# USER SETTINGS
# =========================
subj = "ICE055"

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
    run_feats = sorted(glob.glob(os.path.join(preproc_dir, "Run*.feat")))
    run_feats = [r for r in run_feats if "Exclude" not in os.path.basename(r)]
    if not run_feats:
        raise FileNotFoundError(f"No non-Exclude Run*.feat found in: {preproc_dir}")

    feat_dir = run_feats[0]
    reg_dir = os.path.join(feat_dir, "reg")
    if not os.path.isdir(reg_dir):
        raise FileNotFoundError(f"reg folder not found: {reg_dir}")

    return reg_dir, feat_dir, os.path.basename(feat_dir)


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


def run_std2imgcoord_MNI_to_examplefunc_mm(mni_xyz_mm, reg_dir, out_dir):
    """
    Convert MNI(mm) -> example_func(mm) using:
      std2imgcoord -img example_func -std standard -xfm example_func2standard.mat -mm
    Output: Nx3 mm coords in example_func world space.
    """
    example_func = os.path.join(reg_dir, "example_func.nii.gz")
    standard = os.path.join(reg_dir, "standard.nii.gz")
    xfm = os.path.join(reg_dir, "example_func2standard.mat")

    for p in [example_func, standard, xfm]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required reg file: {p}")

    tmp_in = os.path.join(out_dir, "tmp_mni_coords.txt")
    tmp_out = os.path.join(out_dir, "tmp_examplefunc_mm.txt")
    np.savetxt(tmp_in, mni_xyz_mm, fmt="%.6f")

    cmd = ["std2imgcoord", "-img", example_func, "-std", standard, "-xfm", xfm, "-mm"]

    with open(tmp_in, "r") as fin, open(tmp_out, "w") as fout:
        p = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=subprocess.PIPE, text=True)

    if p.returncode != 0:
        raise RuntimeError(f"std2imgcoord failed:\n{p.stderr}")

    xyz_mm = np.loadtxt(tmp_out)
    if xyz_mm.ndim == 1:
        xyz_mm = xyz_mm.reshape(1, 3)
    if xyz_mm.shape[1] != 3:
        raise ValueError(f"Unexpected output shape from std2imgcoord: {xyz_mm.shape}")

    return xyz_mm, example_func



def affines_close(A, B, tol_mm=2.0):
    # rough check: max absolute difference
    return np.max(np.abs(A - B)) < tol_mm

def mm_to_vox(affine, xyz_mm):
    invA = np.linalg.inv(affine)
    xyz1 = np.c_[xyz_mm, np.ones((xyz_mm.shape[0], 1))]
    ijk = (invA @ xyz1.T).T[:, :3]
    return ijk

def plot_qc(anat_img, pts_ijk, title, out_png, max_points=5000):
    anat = anat_img.get_fdata()
    shp = anat.shape[:3]
    pts = pts_ijk.copy()

    good = np.all(np.isfinite(pts), axis=1)
    pts = pts[good]
    if pts.shape[0] == 0:
        raise ValueError("No valid points to plot.")

    # clip to bounds
    pts[:, 0] = np.clip(pts[:, 0], 0, shp[0]-1)
    pts[:, 1] = np.clip(pts[:, 1], 0, shp[1]-1)
    pts[:, 2] = np.clip(pts[:, 2], 0, shp[2]-1)

    if pts.shape[0] > max_points:
        pts = pts[np.random.choice(pts.shape[0], max_points, replace=False)]

    med = np.median(pts, axis=0).astype(int)
    sx, sy, sz = [int(np.clip(med[i], 0, shp[i]-1)) for i in range(3)]

    dx, dy, dz = anat_img.header.get_zooms()[:3]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Sagittal: x fixed -> show (y,z)
    axes[0].imshow(anat[sx, :, :].T, origin="lower", cmap="gray", aspect=dz/dy)
    axes[0].scatter(pts[:, 1], pts[:, 2], s=10, facecolors='none', edgecolors='r', linewidths=1.2)
    axes[0].set_title(f"Sagittal (x={sx})")
    axes[0].set_xlabel("y (vox)")
    axes[0].set_ylabel("z (vox)")

    # Coronal: y fixed -> show (x,z)
    axes[1].imshow(anat[:, sy, :].T, origin="lower", cmap="gray", aspect=dz/dx)
    axes[1].invert_xaxis()
    axes[1].scatter(pts[:, 0], pts[:, 2], s=10, facecolors='none', edgecolors='r', linewidths=1.2)
    axes[1].set_title(f"Coronal (y={sy})")
    axes[1].set_xlabel("x (vox)")
    axes[1].set_ylabel("z (vox)")

    # Axial: z fixed -> show (x,y)
    axes[2].imshow(anat[:, :, sz].T, origin="lower", cmap="gray", aspect=dy/dx)
    axes[2].invert_xaxis()
    axes[2].scatter(pts[:, 0], pts[:, 1], s=10, facecolors='none', edgecolors='r', linewidths=1.2)
    axes[2].set_title(f"Axial (z={sz})")
    axes[2].set_xlabel("x (vox)")
    axes[2].set_ylabel("y (vox)")

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

# =========================
# MAIN
# =========================
reg_dir, feat_dir, run_used = find_reg_dir(preproc_dir)
print(f"[INFO] Using reg dir: {reg_dir} (from {run_used})")
print(f"[INFO] Using feat dir: {feat_dir}")

# Load Excel MNI coords (mm)
df = pd.read_excel(coords_xlsx)
cx, cy, cz = guess_xyz_cols(df)

mni_xyz = df[[cx, cy, cz]].to_numpy(dtype=float)
good = np.all(np.isfinite(mni_xyz), axis=1)
df_use = df.loc[good].copy()
mni_xyz = mni_xyz[good]
print(f"[INFO] Loaded {len(mni_xyz)} finite MNI coords from Excel.")

# Convert MNI(mm) -> example_func mm coords
func_mm, example_func_nii = run_std2imgcoord_MNI_to_examplefunc_mm(mni_xyz, reg_dir, out_dir)
print("[OK] Converted MNI(mm) -> example_func (mm)")

# ==========================================================
# ADD THIS BLOCK RIGHT HERE
# ==========================================================

# Load functional image
func_img = nib.load(example_func_nii)

# Convert mm → voxel space
func_vox = mm_to_vox(func_img.affine, func_mm)         # float voxel coords
func_vox_round = np.rint(func_vox).astype(int)         # integer voxel indices

# --- build a voxel-space mask with dots at electrode locations ---
mask = np.zeros(func_img.shape[:3], dtype=np.uint8)

ijk = func_vox_round
shp = func_img.shape[:3]

for i, j, k in ijk:
    if 0 <= i < shp[0] and 0 <= j < shp[1] and 0 <= k < shp[2]:
        mask[i, j, k] = 1

mask_nii = nib.Nifti1Image(mask, func_img.affine, func_img.header)
mask_path = os.path.join(out_dir, f"{subj}_electrodes_on_examplefunc_vox_mask.nii.gz")
nib.save(mask_nii, mask_path)
print("[OK] Saved:", mask_path)

# Create 3-panel image (sagittal/coronal/axial)
qc_3panel = os.path.join(out_dir, f"{subj}_QC_ALL_on_examplefunc_3panel.png")
plot_qc(
    func_img,
    func_vox_round,
    title=f"{subj}: ALL electrodes on example_func (3-panel slices)",
    out_png=qc_3panel
)

print("[OK] Saved:", qc_3panel)

# ==========================================================
# THEN CONTINUE WITH YOUR NILEARN OVERLAY
# ==========================================================

qc_func = os.path.join(out_dir, f"{subj}_QC_on_example_func.png")
disp = plotting.plot_epi(example_func_nii,
                         title=f"{subj}: electrodes on example_func (MNI→fMRI mm)")
disp.add_markers(func_mm, marker_size=40, marker_color=(1.0, 0.0, 0.0, 0.35))  # red with alpha=0.35)
disp.savefig(qc_func, dpi=200)
disp.close()
print("[OK] Saved:", qc_func)


# Save CSV
out_csv = os.path.join(out_dir, f"{subj}_MNI_to_examplefunc_mm.csv")
out = df_use.copy()
out["mni_x_mm"] = mni_xyz[:, 0]
out["mni_y_mm"] = mni_xyz[:, 1]
out["mni_z_mm"] = mni_xyz[:, 2]
out["func_x_mm"] = func_mm[:, 0]
out["func_y_mm"] = func_mm[:, 1]
out["func_z_mm"] = func_mm[:, 2]


# --- add voxel coordinates to the output table ---
out["func_i_vox_float"] = func_vox[:, 0]
out["func_j_vox_float"] = func_vox[:, 1]
out["func_k_vox_float"] = func_vox[:, 2]

out["func_i_vox_round"] = np.rint(func_vox[:, 0]).astype(int)
out["func_j_vox_round"] = np.rint(func_vox[:, 1]).astype(int)
out["func_k_vox_round"] = np.rint(func_vox[:, 2]).astype(int)

# bounds check flags (super useful for “are they on the brain?” sanity)
shp = func_img.shape[:3]
out["in_bounds_vox"] = (
    (out["func_i_vox_round"].between(0, shp[0]-1)) &
    (out["func_j_vox_round"].between(0, shp[1]-1)) &
    (out["func_k_vox_round"].between(0, shp[2]-1))
)

# --- SAVE VOXEL COORDS TOO ---
out["func_i_vox_float"] = func_vox[:, 0]
out["func_j_vox_float"] = func_vox[:, 1]
out["func_k_vox_float"] = func_vox[:, 2]

out["func_i_vox_round"] = np.rint(func_vox[:, 0]).astype(int)
out["func_j_vox_round"] = np.rint(func_vox[:, 1]).astype(int)
out["func_k_vox_round"] = np.rint(func_vox[:, 2]).astype(int)

# optional: quick bounds flag
shp = func_img.shape[:3]
out["in_bounds_vox"] = (
    (out["func_i_vox_round"].between(0, shp[0]-1)) &
    (out["func_j_vox_round"].between(0, shp[1]-1)) &
    (out["func_k_vox_round"].between(0, shp[2]-1))
)


out.to_csv(out_csv, index=False)
print("[OK] Saved:", out_csv)

print("[DONE]")
