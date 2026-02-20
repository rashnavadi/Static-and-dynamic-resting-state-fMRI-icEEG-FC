import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nilearn import plotting
import re
import numpy.linalg as npl

# -------------------------
# PATHS
# -------------------------
seed_csv = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/python_inputs/ICE060/ICE060_Run1_seedNames.csv"
coords_csv = "/Volumes/MIND/ICE/scripts/Kristina_paper/ICE060_native_space_coords.csv"
anat_native = "/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/ICE060/4_MRI/1_Anatomicals/3D_Anat/ICE060_3Danat_brain.nii.gz"

out_png_ortho = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/ICE060_electrodes_native_ortho.png"
out_png_yz    = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/ICE060_electrodes_native_yz.png"

# -------------------------
# Helper functions
# -------------------------
def canonical_contact(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip().replace(" ", "").replace("_", "")
    s = s.upper()
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s

def parse_bipolar(seed: str):
    s = str(seed).strip().upper()
    if "-" not in s:
        return None, None
    a, b = s.split("-", 1)
    return canonical_contact(a), canonical_contact(b)

def world_to_vox(xyz_mm, inv_aff):
    xyz_mm = np.asarray(xyz_mm, float)
    ones = np.ones((xyz_mm.shape[0], 1))
    xyz1 = np.hstack([xyz_mm, ones])
    ijk1 = xyz1 @ inv_aff.T
    return ijk1[:, :3]

def frac_in_bounds(ijk, shape):
    ijk = np.asarray(ijk)
    ok = (
        (ijk[:, 0] >= 0) & (ijk[:, 0] < shape[0]) &
        (ijk[:, 1] >= 0) & (ijk[:, 1] < shape[1]) &
        (ijk[:, 2] >= 0) & (ijk[:, 2] < shape[2])
    )
    return ok.mean()

def make_markers_hollow(display, edgecolor="red", alpha=0.7, lw=1.2):
    """
    Convert the *last added* marker PathCollections to hollow markers.
    Works on nilearn displays (plot_anat, plot_glass_brain, etc.).
    """
    for cutax in display.axes.values():      # CutAxes
        ax = cutax.ax                        # matplotlib Axes
        if len(ax.collections) > 0:
            coll = ax.collections[-1]        # last added = our markers
            coll.set_facecolor("none")       # hollow
            coll.set_edgecolor(edgecolor)
            coll.set_linewidth(lw)
            coll.set_alpha(alpha)

# -------------------------
# 1) Load seeds
# -------------------------
print(">>> Loading seedNames CSV")
seed_df = pd.read_csv(seed_csv)
seed_names = [str(s).strip().upper() for s in seed_df.iloc[:, 0].tolist() if str(s).strip() != ""]
print(f"[INFO] seedNames N={len(seed_names)}")
print("[INFO] first 10:", seed_names[:10])

# -------------------------
# 2) Load coords and build monopolar contact names
# -------------------------
print(">>> Loading native coords CSV")
df = pd.read_csv(coords_csv)
print(f"[INFO] coords CSV shape: {df.shape}")
print("[INFO] coords CSV columns:", list(df.columns))

shaft_col = "ElectrodeName_labConvention_"
num_col   = "ContactNumber"
xcol, ycol, zcol = "x", "y", "z"

df = df.dropna(subset=[shaft_col, num_col, xcol, ycol, zcol]).copy()
df["_mono"] = df[shaft_col].astype(str).apply(canonical_contact) + df[num_col].astype(int).astype(str)

# contact -> xyz (whatever units are in CSV)
contact_to_xyz = {}
for _, row in df.iterrows():
    contact_to_xyz[row["_mono"]] = np.array([row[xcol], row[ycol], row[zcol]], float)

# -------------------------
# 3) Make bipolar midpoints
# -------------------------
coords_raw = []
matched_seeds = []
missing = []

for bn in seed_names:
    a, b = parse_bipolar(bn)
    if a is None:
        missing.append((bn, "not bipolar"))
        continue
    if a not in contact_to_xyz or b not in contact_to_xyz:
        missing.append((bn, f"missing {a if a not in contact_to_xyz else ''} {b if b not in contact_to_xyz else ''}".strip()))
        continue
    coords_raw.append(0.5 * (contact_to_xyz[a] + contact_to_xyz[b]))
    matched_seeds.append(bn)

coords_raw = np.asarray(coords_raw, float)
print(f"[INFO] Matched bipolar seeds: {len(matched_seeds)} / {len(seed_names)}")
if missing:
    print("[WARN] Missing (first 10):", missing[:10])

if len(coords_raw) == 0:
    raise RuntimeError("No bipolar midpoints were created. Check your naming columns.")

# -------------------------
# 4) Check if coords_raw behave like WORLD(mm) in this T1 space
#    (Your prints suggest YES for ICE060)
# -------------------------
anat_img = nib.load(anat_native)
aff = anat_img.affine
inv_aff = npl.inv(aff)
shape = anat_img.shape[:3]

ijk_if_world = world_to_vox(coords_raw, inv_aff)
frac_world = frac_in_bounds(ijk_if_world, shape)

print("\n[DEBUG] Anatomy shape:", shape)
print("[DEBUG] Anatomy affine:\n", aff)
print("[DEBUG] Treat coords_raw as WORLD(mm) -> frac in bounds:", f"{frac_world:.2f}")
print("[DEBUG] ijk_if_world min/max:", ijk_if_world.min(axis=0), ijk_if_world.max(axis=0))

# For ICE060 you got ~1.00, so use as world/mm
coords_world = coords_raw.copy()
print("[INFO] Using coords as WORLD(mm) (native).")
print("[INFO] coords_world min/max:", coords_world.min(axis=0), coords_world.max(axis=0))

# -------------------------
# 5) ORTHO plot (FreeView-like)
# -------------------------
cut_xyz = np.median(coords_world, axis=0).tolist()
print("[INFO] cut_xyz (median):", [round(v, 3) for v in cut_xyz])

marker_size = 18     # smaller dots
alpha = 0.7          # transparency
lw = 1.2             # outline thickness

print(">>> Plotting ORTHO + saving:", out_png_ortho)
disp = plotting.plot_anat(
    anat_native,
    display_mode="ortho",
    cut_coords=cut_xyz,
    title="ICE060 bipolar seeds (native ORTHO)",
    annotate=True,
    draw_cross=False,
)
disp.add_markers(coords_world, marker_color="red", marker_size=marker_size)
make_markers_hollow(disp, edgecolor="red", alpha=alpha, lw=lw)
disp.savefig(out_png_ortho)
disp.close()
print("[OK] Saved:", out_png_ortho)

# -------------------------
# 6) Y / Z plot (1x2) using matplotlib axes
# -------------------------
cut_y = float(np.median(coords_world[:, 1]))
cut_z = float(np.median(coords_world[:, 2]))
print("[INFO] cut_y, cut_z:", round(cut_y, 3), round(cut_z, 3))

print(">>> Plotting Y/Z + saving:", out_png_yz)
fig = plt.figure(figsize=(11, 5))

ax1 = fig.add_subplot(1, 2, 1)
d1 = plotting.plot_anat(
    anat_native, display_mode="y", cut_coords=[cut_y], axes=ax1,
    title="Coronal (y)", annotate=True, draw_cross=False
)
d1.add_markers(coords_world, marker_color="red", marker_size=marker_size)
make_markers_hollow(d1, edgecolor="red", alpha=alpha, lw=lw)

ax2 = fig.add_subplot(1, 2, 2)
d2 = plotting.plot_anat(
    anat_native, display_mode="z", cut_coords=[cut_z], axes=ax2,
    title="Axial (z)", annotate=True, draw_cross=False
)
d2.add_markers(coords_world, marker_color="red", marker_size=marker_size)
make_markers_hollow(d2, edgecolor="red", alpha=alpha, lw=lw)

plt.tight_layout()
plt.savefig(out_png_yz, dpi=300, bbox_inches="tight")
plt.close(fig)
print("[OK] Saved:", out_png_yz)

print(">>> DONE")
