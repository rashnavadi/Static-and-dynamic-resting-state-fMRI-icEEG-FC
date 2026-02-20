import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import re
from nilearn import plotting


# -------------------------
# PATHS (EDIT IF NEEDED)
# -------------------------
seed_csv = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/python_inputs/ICE060/ICE060_Run1_seedNames.csv"

# Labeled FC CSV exported from MATLAB (with row+col names)
fc_csv_labeled = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/ICE060/ICE060_Run1_FC_static_labeled.csv"

coords_csv = "/Volumes/MIND/ICE/scripts/Kristina_paper/ICE060_native_space_coords.csv"
anat_native = "/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/ICE060/4_MRI/1_Anatomicals/3D_Anat/ICE060_3Danat_brain.nii.gz"

out_png = "/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/ICE060_top20_edges_connectome.png"

TOP_K = 20


# -------------------------
# Helpers
# -------------------------
def canonical(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip().replace(" ", "").replace("_", "")
    s = s.upper()
    s = re.sub(r"[^A-Z0-9\-]", "", s)
    return s

def canonical_contact(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip().replace(" ", "").replace("_", "")
    s = s.upper()
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s

def parse_bipolar(seed: str):
    s = canonical(seed)
    if "-" not in s:
        return None, None
    a, b = s.split("-", 1)
    return canonical_contact(a), canonical_contact(b)

def top_k_edges_by_abs(A, k=20):
    """Return list of (i,j,val) for top-k abs edges in upper triangle (excluding diag)."""
    A = np.asarray(A, float)
    n = A.shape[0]
    iu = np.triu_indices(n, k=1)
    vals = A[iu]
    order = np.argsort(np.abs(vals))[::-1]
    order = order[: min(k, len(order))]
    edges = []
    for idx in order:
        i = iu[0][idx]
        j = iu[1][idx]
        edges.append((i, j, A[i, j]))
    return edges

def family_label(seed_name: str) -> str:
    # 'DLA1-DLA2' -> 'DLA' ; 'DLAH1-DLAH2' -> 'DLAH'
    s = canonical(seed_name).replace("-", "")
    s = re.sub(r"\d+", "", s)  # remove digits
    return s

def select_farthest_edges(node_coords, seed_order, anat_native, k=15,
                          force_cross_hemi=True, force_diff_family=True):
    """
    Returns list of (i, j, dist_mm) for top-k farthest node pairs.
    Uses voxel i-index (after world->vox) to determine hemisphere robustly.
    """
    coords = np.asarray(node_coords, float)
    N = coords.shape[0]

    # hemisphere via voxel i-index
    img = nib.load(anat_native)
    inv_aff = np.linalg.inv(img.affine)
    shape = img.shape[:3]
    mid_i = shape[0] / 2.0

    ones = np.ones((N, 1))
    ijk = (np.hstack([coords, ones]) @ inv_aff.T)[:, :3]
    hemi = np.where(ijk[:, 0] < mid_i, "L", "R")

    pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            if force_cross_hemi and (hemi[i] == hemi[j]):
                continue
            if force_diff_family:
                if family_label(seed_order[i]) == family_label(seed_order[j]):
                    continue
            dist = float(np.linalg.norm(coords[i] - coords[j]))  # mm
            pairs.append((dist, i, j))

    pairs.sort(reverse=True, key=lambda x: x[0])
    pairs = pairs[: min(k, len(pairs))]

    return pairs, hemi

# -------------------------
# 1) Load seed order
# -------------------------
print(">>> Loading seed order")
seed_df = pd.read_csv(seed_csv)
seed_order = [canonical(x) for x in seed_df.iloc[:, 0].tolist() if str(x).strip() != ""]
print(f"[INFO] seedNames N={len(seed_order)}")
print("[INFO] first 10:", seed_order[:10])


# -------------------------
# 2) Load FC matrix (labeled) and reorder to match seed_order
# -------------------------
print(">>> Loading FC CSV (labeled)")
df_fc = pd.read_csv(fc_csv_labeled, index_col=0)
df_fc.index = [canonical(x) for x in df_fc.index]
df_fc.columns = [canonical(x) for x in df_fc.columns]

# Reindex to match seed_order (this guarantees FC matches your node coords)
missing_rows = [s for s in seed_order if s not in df_fc.index]
missing_cols = [s for s in seed_order if s not in df_fc.columns]
if missing_rows or missing_cols:
    raise ValueError(
        f"FC CSV labels do not match seed order.\n"
        f"Missing rows: {missing_rows[:10]}\n"
        f"Missing cols: {missing_cols[:10]}"
    )

df_fc = df_fc.loc[seed_order, seed_order]
FC = df_fc.values.astype(float)
print("[INFO] FC shape:", FC.shape)


# -------------------------
# 3) Load native coords and build monopolar -> xyz
#    monopolar = ElectrodeName_labConvention_ + ContactNumber  (e.g., dLaH + 3 -> DLAH3)
# -------------------------
print(">>> Loading coords CSV and building monopolar map")
df = pd.read_csv(coords_csv)

shaft_col = "ElectrodeName_labConvention_"
num_col   = "ContactNumber"
xcol, ycol, zcol = "x", "y", "z"

for col in [shaft_col, num_col, xcol, ycol, zcol]:
    if col not in df.columns:
        raise ValueError(f"Missing column '{col}' in coords CSV.")

df = df.dropna(subset=[xcol, ycol, zcol]).copy()
df["_mono"] = df[shaft_col].astype(str).apply(canonical_contact) + df[num_col].astype(int).astype(str)

contact_to_xyz = {
    row["_mono"]: np.array([row[xcol], row[ycol], row[zcol]], float)
    for _, row in df.iterrows()
}

# Build bipolar midpoints in WORLD(mm) native space (your CSV appears to be world coords)
node_coords = []
kept = []
for seed in seed_order:
    a, b = parse_bipolar(seed)
    if a not in contact_to_xyz or b not in contact_to_xyz:
        raise ValueError(f"Missing monopolar contacts for seed {seed}: {a}, {b}")
    mid = 0.5 * (contact_to_xyz[a] + contact_to_xyz[b])
    node_coords.append(mid)
    kept.append(seed)

node_coords = np.asarray(node_coords, float)
print("[INFO] node_coords shape:", node_coords.shape)
print("[INFO] node_coords min/max:", node_coords.min(axis=0), node_coords.max(axis=0))


# -------------------------
# 4) Select TOP_K strongest edges (abs) and build sparse adjacency
# -------------------------
# -------------------------
# 4) Select farthest electrode pairs (NOT strongest FC) and build adjacency
# -------------------------
TOP_K = 15  # keep 10–20 for richer figure

pairs, hemi = select_farthest_edges(
    node_coords=node_coords,
    seed_order=seed_order,
    anat_native=anat_native,
    k=TOP_K,
    force_cross_hemi=True,     # show long-range edges across hemispheres
    force_diff_family=True     # avoid edges within same shaft family
)

A_top = np.zeros_like(FC)

# Option A: draw edges with equal weight (simple demo)
for dist, i, j in pairs:
    A_top[i, j] = 1.0
    A_top[j, i] = 1.0

# Option B (if you prefer edge thickness reflects FC at those far distances):
# for dist, i, j in pairs:
#     A_top[i, j] = FC[i, j]
#     A_top[j, i] = FC[i, j]

print(f"[INFO] Selected top {len(pairs)} edges by DISTANCE (mm):")
for dist, i, j in pairs[:10]:
    print(f"   {seed_order[i]} ({hemi[i]}) -- {seed_order[j]} ({hemi[j]})   dist={dist:.1f} mm   FC={FC[i,j]:.4f}")

# -------------------------
# 5) Plot connectome (glass brain)
# -------------------------
print(">>> Plotting connectome + saving:", out_png)

disp = plotting.plot_connectome(
    A_top,
    node_coords,
    node_size=12,            # smaller nodes
    display_mode="ortho",
    title=f"ICE060 farthest {len(pairs)} electrode pairs"
)

disp.savefig(out_png)
disp.close()

print("[OK] Saved:", out_png)
print(">>> DONE")
