#!/usr/bin/env python3
"""
ICE041 native-space electrode localization on anatomical MRI.

Inputs:
  - anat_nii: ICE041_3Danat_brain.nii.gz (native space)
  - coords_xlsx: ICE041_native_space_coords.xlsx (native space coords)

Outputs:
  - electrodes_localized_native.csv
  - QC PNGs with electrode overlays

Optional:
  - label_nii: parcellation/segmentation NIfTI in SAME native space as anat
  - lut_csv: label lookup table (label integer -> name)
"""

import os
import re
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import nibabel as nib


# -----------------------
# User paths (EDIT IF NEEDED)
# -----------------------
anat_nii = "/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/ICE060/4_MRI/1_Anatomicals/3D_Anat/ICE060_3Danat_brain.nii.gz"
coords_xlsx = "/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/ICE060_native_space_coords.xlsx"

# Optional: native-space label map + LUT
label_nii = None  # e.g., "/Volumes/.../ICE041_aparc+aseg_native.nii.gz"
lut_csv   = None  # e.g., "/Volumes/.../FreeSurferColorLUT.csv" or your own CSV with columns: id,name

out_dir = "/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/ICE060_outputs"
os.makedirs(out_dir, exist_ok=True)

# -----------------------
# Helpers
# -----------------------
def guess_coord_columns(df: pd.DataFrame):
    """
    Tries to find x/y/z columns (in mm) automatically.
    Accepts common patterns: x y z, X Y Z, coord_x, mni_x, native_x, etc.
    """
    cols = list(df.columns)

    def find_one(axis):
        pats = [
            rf"^{axis}$",
            rf"^{axis}\s*\(mm\)$",
            rf".*{axis}.*mm.*",
            rf".*coord.*{axis}.*",
            rf".*native.*{axis}.*",
            rf".*mni.*{axis}.*",
        ]
        for p in pats:
            for c in cols:
                if re.match(p, str(c).strip(), flags=re.IGNORECASE):
                    return c
        return None

    cx, cy, cz = find_one("x"), find_one("y"), find_one("z")

    # fallback: look for 3 numeric columns that appear together
    if cx is None or cy is None or cz is None:
        num_cols = [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]
        # try to find consecutive triplets that look like coordinates
        for i in range(len(num_cols) - 2):
            a, b, c = num_cols[i:i+3]
            # heuristic: coords usually have both positive/negative values in mm space
            if (df[a].abs().median() > 1 and df[a].abs().median() < 300 and
                df[b].abs().median() > 1 and df[b].abs().median() < 300 and
                df[c].abs().median() > 1 and df[c].abs().median() < 300):
                cx, cy, cz = a, b, c
                break

    if cx is None or cy is None or cz is None:
        raise ValueError(
            "Could not auto-detect coordinate columns. "
            "Please rename your coord columns to x,y,z (in mm) or edit guess_coord_columns()."
        )

    return cx, cy, cz

def world_to_voxel(affine, xyz_mm):
    """Convert world mm (x,y,z) to voxel (i,j,k) using NIfTI affine."""
    invA = np.linalg.inv(affine)
    xyz1 = np.c_[xyz_mm, np.ones((xyz_mm.shape[0], 1))]
    ijk  = (invA @ xyz1.T).T[:, :3]
    return ijk

def voxel_to_world(affine, ijk):
    """
    Convert voxel coords to world (mm).
    Supports:
      - ijk shape (3,)   -> returns (3,)
      - ijk shape (N,3)  -> returns (N,3)
    """
    ijk = np.asarray(ijk)
    if ijk.ndim == 1:
        ijk1 = np.r_[ijk, 1.0]                      # (4,)
        return (affine @ ijk1)[:3]                  # (3,)
    elif ijk.ndim == 2 and ijk.shape[1] == 3:
        ijk1 = np.c_[ijk, np.ones((ijk.shape[0],1))]  # (N,4)
        return (affine @ ijk1.T).T[:, :3]             # (N,3)
    else:
        raise ValueError(f"ijk must be shape (3,) or (N,3), got {ijk.shape}")


def sample_label(label_img, ijk_round):
    """Sample integer label volume at voxel coords (rounded/clipped)."""
    data = label_img.get_fdata()
    shp = data.shape
    ijk = ijk_round.astype(int)
    ijk[:, 0] = np.clip(ijk[:, 0], 0, shp[0]-1)
    ijk[:, 1] = np.clip(ijk[:, 1], 0, shp[1]-1)
    ijk[:, 2] = np.clip(ijk[:, 2], 0, shp[2]-1)
    vals = data[ijk[:, 0], ijk[:, 1], ijk[:, 2]].astype(int)
    return vals

def load_lut(lut_path):
    """
    Load LUT with columns: id,name
    Accepts CSV/TSV. If FreeSurferColorLUT.txt, you’ll want to convert to CSV first.
    """
    lut = pd.read_csv(lut_path)
    if "id" not in lut.columns or "name" not in lut.columns:
        raise ValueError("LUT must have columns: id,name")
    return dict(zip(lut["id"].astype(int), lut["name"].astype(str)))

def plot_qc(anat_img, points_ijk, title, out_png, max_points=400):
    anat = anat_img.get_fdata()
    dx, dy, dz = anat_img.header.get_zooms()[:3]  # (dx, dy, dz) in mm
    shp = anat.shape

    pts = points_ijk.copy()
    if pts.shape[0] > max_points:
        pts = pts[np.random.choice(pts.shape[0], max_points, replace=False)]

    med = np.median(pts, axis=0).astype(int)
    sx, sy, sz = [int(np.clip(med[i], 0, shp[i]-1)) for i in range(3)]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    mid_x_vox = shp[0] // 2
    axes[1].axvline(mid_x_vox, linewidth=1)
    axes[2].axvline(mid_x_vox, linewidth=1)
    # Sagittal: x fixed -> image is (y,z), but we plot .T so rows=z, cols=y
    axes[0].imshow(anat[sx, :, :].T, origin="lower", cmap="gray", aspect=dz/dy)
    axes[0].scatter(pts[:, 1], pts[:, 2], s=10)
    axes[0].set_title(f"Sagittal (x={sx})")
    axes[0].set_xlabel("y (vox)")
    axes[0].set_ylabel("z (vox)")

    # Coronal: y fixed -> image is (x,z), .T => rows=z, cols=x
    axes[1].imshow(anat[:, sy, :].T, origin="lower", cmap="gray", aspect=dz/dx)
    axes[1].invert_xaxis()  # coronal: x axis
    axes[1].scatter(pts[:, 0], pts[:, 2], s=10)
    axes[1].set_title(f"Coronal (y={sy})")
    axes[1].set_xlabel("x (vox)")
    axes[1].set_ylabel("z (vox)")

    # Axial: z fixed -> image is (x,y), .T => rows=y, cols=x
    axes[2].imshow(anat[:, :, sz].T, origin="lower", cmap="gray", aspect=dy/dx)
    axes[2].invert_xaxis()  # axial: x axis
    axes[2].scatter(pts[:, 0], pts[:, 1], s=10)
    axes[2].set_title(f"Axial (z={sz})")
    axes[2].set_xlabel("x (vox)")
    axes[2].set_ylabel("y (vox)")

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)

# -----------------------
# Main
# -----------------------
print("Loading anatomy:", anat_nii)
anat_img = nib.load(anat_nii)
print("Anat orientation (aff2axcodes):", nib.aff2axcodes(anat_img.affine))
print("Voxel sizes (mm):", anat_img.header.get_zooms()[:3])
print("Affine:\n", anat_img.affine)

shp = anat_img.shape[:3]
A   = anat_img.affine

jmid, kmid = shp[1]//2, shp[2]//2

x0   = voxel_to_world(A, [0,        jmid, kmid])[0]
xend = voxel_to_world(A, [shp[0]-1, jmid, kmid])[0]

x_min = min(x0, xend)
x_max = max(x0, xend)
x_mid = 0.5*(x_min + x_max)

print("World-x range:", x_min, "to", x_max)
print("Estimated mid-sagittal x_mid:", x_mid)

print("Loading Excel:", coords_xlsx)
df = pd.read_excel(coords_xlsx)

# --- required columns from your message ---
# column 6: ElectrodeName_labConvention_
# column 3: number of contacts
# plus IncludedInIEEG_fMRIStudy
# We’ll locate by name rather than hard-coded index, but keep your intent.

name_col_candidates = [c for c in df.columns if "ElectrodeName_labConvention" in str(c)]
if len(name_col_candidates) == 0:
    # fallback: any column containing "ElectrodeName" and "lab"
    name_col_candidates = [c for c in df.columns if ("electrode" in str(c).lower() and "name" in str(c).lower())]
if len(name_col_candidates) == 0:
    raise ValueError("Could not find ElectrodeName_labConvention_ column in Excel.")
name_col = name_col_candidates[0]

included_candidates = [c for c in df.columns if "IncludedInIEEG_fMRIStudy" in str(c)]
if len(included_candidates) == 0:
    raise ValueError("Could not find IncludedInIEEG_fMRIStudy column in Excel.")
included_col = included_candidates[0]

# contact count (column 3). We try to find a sensible column name:
contact_candidates = [c for c in df.columns if ("contact" in str(c).lower() and "num" in str(c).lower())]
if len(contact_candidates) == 0:
    # fallback: any column with "contact" and numeric
    contact_candidates = [c for c in df.columns if ("contact" in str(c).lower() and pd.api.types.is_numeric_dtype(df[c]))]
if len(contact_candidates) == 0:
    # last resort: use 3rd column index (your note)
    contact_col = df.columns[2]
else:
    contact_col = contact_candidates[0]

# coordinate columns
cx, cy, cz = guess_coord_columns(df)

# filter included
df_use = df[df[included_col] == 1].copy()
df_use = df_use.dropna(subset=[cx, cy, cz])
print(f"Included rows: {len(df_use)}")

# Extract xyz (assumed mm in native/world space of the NIfTI)
xyz = df_use[[cx, cy, cz]].to_numpy(dtype=float)

# Convert to voxel coords
ijk = world_to_voxel(A, xyz)
ijk_round = np.rint(ijk).astype(int)

# Convert back to world (for sanity check; should be close)
xyz_back = voxel_to_world(A, ijk_round)

# Optional label sampling
region_id = None
region_name = None
if label_nii is not None:
    print("Loading label map:", label_nii)
    lab_img = nib.load(label_nii)
    # sanity: must match grid/affine reasonably
    if lab_img.shape[:3] != anat_img.shape[:3]:
        print("[WARN] label_nii shape != anat shape. They must be in the SAME native space/grid for correct lookup.")
    region_id = sample_label(lab_img, ijk_round)
    if lut_csv is not None:
        lut = load_lut(lut_csv)
        region_name = [lut.get(int(i), f"Unknown_{int(i)}") for i in region_id]

# Build output table
out = df_use.copy()
out["coord_space"] = "native_mm_assumed"
out["x_mm"] = xyz[:, 0]
out["y_mm"] = xyz[:, 1]
out["z_mm"] = xyz[:, 2]
out["i_vox"] = ijk_round[:, 0]
out["j_vox"] = ijk_round[:, 1]
out["k_vox"] = ijk_round[:, 2]
out["x_mm_from_vox"] = xyz_back[:, 0]
out["y_mm_from_vox"] = xyz_back[:, 1]
out["z_mm_from_vox"] = xyz_back[:, 2]

if region_id is not None:
    out["region_id"] = region_id
if region_name is not None:
    out["region_name"] = region_name

for prefix in ["dLA", "dRA"]:
    m = out[name_col].astype(str).str.startswith(prefix)
    if m.any():
        print(prefix, "n=", m.sum(),
              "mean i_vox=", out.loc[m, "i_vox"].mean(),
              "mean x_mm=", out.loc[m, "x_mm"].mean())

for prefix in ["dLA", "dRA"]:
    m = out[name_col].astype(str).str.startswith(prefix)
    if m.any():
        xs = out.loc[m, "x_mm"].to_numpy()
        left_count  = int(np.sum(xs > x_mid))
        right_count = int(np.sum(xs <= x_mid))
        print(prefix, f"LEFT_of_mid={left_count}, RIGHT_of_mid={right_count}")

# Save CSV
csv_path = os.path.join(out_dir, "electrodes_localized_native.csv")
out.to_csv(csv_path, index=False)
print("Saved:", csv_path)

# QC plots (all included points)
qc_png = os.path.join(out_dir, "QC_electrodes_on_anat.png")
plot_qc(anat_img, ijk_round, "ICE060 electrodes (IncludedInStudy==1) on native anatomy", qc_png)
print("Saved:", qc_png)

# QC plots grouped by ElectrodeName_labConvention_ (shaft/lead)
# (helps if each row is a contact and name indicates shaft)
for shaft, g in out.groupby(name_col):
    pts = g[["i_vox","j_vox","k_vox"]].to_numpy(int)
    if pts.shape[0] < 2:
        continue
    safe = re.sub(r"[^A-Za-z0-9_\-]+", "_", str(shaft))
    png = os.path.join(out_dir, f"QC_{safe}.png")
    plot_qc(anat_img, pts, f"{shaft} (n={pts.shape[0]})", png, max_points=999999)

print("Done.")
