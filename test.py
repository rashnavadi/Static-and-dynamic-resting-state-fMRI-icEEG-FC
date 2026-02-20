#!/usr/bin/env python3
import os
import re
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt


# =========================
# USER PATHS (EDIT)
# =========================
subj = "ICE013"

anat_native_nii = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}/4_MRI/1_Anatomicals/3D_Anat/{subj}_3Danat_brain.nii.gz"

reg_dir = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}/4_MRI/2_Functionals/2_PreProcessing/Run1a.feat/reg"

# MNI coord Excel (contains x,y,z in MNI mm)
coords_xlsx = f"/Volumes/MIND/ICE/original_ICE/MNI_coordinates/{subj}_channel_info.xlsx"

out_dir = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_outputs_py"
os.makedirs(out_dir, exist_ok=True)


# =========================
# HELPERS
# =========================
def guess_coord_columns(df: pd.DataFrame):
    cols = list(df.columns)

    def norm(s):
        return re.sub(r"[^a-z0-9]+", "", str(s).strip().lower())

    nmap = {c: norm(c) for c in cols}

    # prefer explicit mni_* if present
    candidates = {
        "x": ["mnix", "mni_x", "x", "xmm", "coordx", "posx"],
        "y": ["mniy", "mni_y", "y", "ymm", "coordy", "posy"],
        "z": ["mniz", "mni_z", "z", "zmm", "coordz", "posz"],
    }

    def pick(axis):
        for c in cols:
            if nmap[c] in [norm(v) for v in candidates[axis]]:
                return c
        # fallback: contains axis and mni
        for c in cols:
            s = nmap[c]
            if ("mni" in s) and (axis in s):
                return c
        # last resort: exact axis
        for c in cols:
            if nmap[c] == axis:
                return c
        return None

    cx, cy, cz = pick("x"), pick("y"), pick("z")
    if cx is None or cy is None or cz is None:
        raise ValueError(f"Could not detect x/y/z columns. Columns: {cols}")
    return cx, cy, cz


def load_fsl_flirt_mat(mat_path: str) -> np.ndarray:
    """Load a 4x4 FLIRT matrix saved as text (usually 4 lines of 4 floats)."""
    M = np.loadtxt(mat_path)
    if M.shape != (4, 4):
        raise ValueError(f"Expected 4x4 matrix in {mat_path}, got {M.shape}")
    return M


def mm_to_vox(affine: np.ndarray, xyz_mm: np.ndarray) -> np.ndarray:
    """World mm -> voxel (float) using NIfTI affine."""
    invA = np.linalg.inv(affine)
    xyz1 = np.c_[xyz_mm, np.ones((xyz_mm.shape[0], 1))]
    ijk = (invA @ xyz1.T).T[:, :3]
    return ijk


def vox_to_mm(affine: np.ndarray, ijk: np.ndarray) -> np.ndarray:
    """Voxel (float) -> world mm using NIfTI affine."""
    ijk1 = np.c_[ijk, np.ones((ijk.shape[0], 1))]
    xyz = (affine @ ijk1.T).T[:, :3]
    return xyz


def apply_flirt_vox2vox(M_vox2vox: np.ndarray, ijk_src: np.ndarray) -> np.ndarray:
    """
    Apply FLIRT matrix in voxel space:
      ijk_tgt_h = M * ijk_src_h
    """
    ijk1 = np.c_[ijk_src, np.ones((ijk_src.shape[0], 1))]
    out = (M_vox2vox @ ijk1.T).T
    return out[:, :3]


def plot_qc(anat_img, pts_ijk, title, out_png, max_points=5000):
    anat = anat_img.get_fdata()
    shp = anat.shape[:3]
    pts = pts_ijk.copy()

    # clip + remove nan/inf
    good = np.all(np.isfinite(pts), axis=1)
    pts = pts[good]

    if pts.shape[0] == 0:
        raise ValueError("No valid points to plot.")

    # clip to bounds for plotting
    pts[:, 0] = np.clip(pts[:, 0], 0, shp[0]-1)
    pts[:, 1] = np.clip(pts[:, 1], 0, shp[1]-1)
    pts[:, 2] = np.clip(pts[:, 2], 0, shp[2]-1)

    if pts.shape[0] > max_points:
        pts = pts[np.random.choice(pts.shape[0], max_points, replace=False)]

    med = np.median(pts, axis=0).astype(int)
    sx, sy, sz = [int(np.clip(med[i], 0, shp[i]-1)) for i in range(3)]

    dx, dy, dz = anat_img.header.get_zooms()[:3]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Sagittal x fixed: show (y,z)
    axes[0].imshow(anat[sx, :, :].T, origin="lower", cmap="gray", aspect=dz/dy)
    axes[0].scatter(pts[:, 1], pts[:, 2], s=10)
    axes[0].set_title(f"Sagittal (x={sx})")
    axes[0].set_xlabel("y (vox)")
    axes[0].set_ylabel("z (vox)")

    # Coronal y fixed: show (x,z)
    axes[1].imshow(anat[:, sy, :].T, origin="lower", cmap="gray", aspect=dz/dx)
    axes[1].invert_xaxis()
    axes[1].scatter(pts[:, 0], pts[:, 2], s=10)
    axes[1].set_title(f"Coronal (y={sy})")
    axes[1].set_xlabel("x (vox)")
    axes[1].set_ylabel("z (vox)")

    # Axial z fixed: show (x,y)
    axes[2].imshow(anat[:, :, sz].T, origin="lower", cmap="gray", aspect=dy/dx)
    axes[2].invert_xaxis()
    axes[2].scatter(pts[:, 0], pts[:, 1], s=10)
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
# reg-space images + matrix
standard_nii = os.path.join(reg_dir, "standard.nii.gz")
highres_nii  = os.path.join(reg_dir, "highres.nii.gz")
M_std2hi_mat = os.path.join(reg_dir, "standard2highres.mat")

for p in [standard_nii, highres_nii, M_std2hi_mat, anat_native_nii, coords_xlsx]:
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing: {p}")

print("[INFO] Loading images/matrix...")
std_img = nib.load(standard_nii)
hi_img  = nib.load(highres_nii)
anat_img = nib.load(anat_native_nii)
M_std2hi = load_fsl_flirt_mat(M_std2hi_mat)

print("[INFO] standard aff2axcodes:", nib.aff2axcodes(std_img.affine))
print("[INFO] highres   aff2axcodes:", nib.aff2axcodes(hi_img.affine))
print("[INFO] anat      aff2axcodes:", nib.aff2axcodes(anat_img.affine))

print("[INFO] Reading Excel:", coords_xlsx)
df = pd.read_excel(coords_xlsx)
cx, cy, cz = guess_coord_columns(df)

xyz_mni_mm = df[[cx, cy, cz]].to_numpy(dtype=float)

# drop rows with NaNs in coords
good = np.all(np.isfinite(xyz_mni_mm), axis=1)
df_use = df.loc[good].copy()
xyz_mni_mm = xyz_mni_mm[good]

print(f"[INFO] Using {xyz_mni_mm.shape[0]} coordinate rows (finite xyz).")

# --- Step 1: MNI mm -> standard vox
ijk_std = mm_to_vox(std_img.affine, xyz_mni_mm)

# --- Step 2: standard vox -> highres vox (apply FLIRT std2highres)
ijk_hi = apply_flirt_vox2vox(M_std2hi, ijk_std)

# --- Step 3: highres vox -> highres mm
xyz_hi_mm = vox_to_mm(hi_img.affine, ijk_hi)

# --- Step 4: highres mm -> native anat vox (overlay target)
ijk_anat = mm_to_vox(anat_img.affine, xyz_hi_mm)
ijk_anat_round = np.rint(ijk_anat).astype(int)

# Save output table
out = df_use.copy()
out["mni_x_mm"] = xyz_mni_mm[:, 0]
out["mni_y_mm"] = xyz_mni_mm[:, 1]
out["mni_z_mm"] = xyz_mni_mm[:, 2]

out["highres_x_mm"] = xyz_hi_mm[:, 0]
out["highres_y_mm"] = xyz_hi_mm[:, 1]
out["highres_z_mm"] = xyz_hi_mm[:, 2]

out["anat_i_vox"] = ijk_anat_round[:, 0]
out["anat_j_vox"] = ijk_anat_round[:, 1]
out["anat_k_vox"] = ijk_anat_round[:, 2]

csv_out = os.path.join(out_dir, f"{subj}_coords_MNI_to_nativeAnat.csv")
out.to_csv(csv_out, index=False)
print("[OK] Saved:", csv_out)

# QC plot on native anat
qc_png = os.path.join(out_dir, f"{subj}_QC_MNIcoords_on_nativeAnat.png")
plot_qc(anat_img, ijk_anat_round.astype(float),
        title=f"{subj}: MNI coords transformed to native anat (std2highres.mat)",
        out_png=qc_png)
print("[OK] Saved:", qc_png)

print("[DONE]")
