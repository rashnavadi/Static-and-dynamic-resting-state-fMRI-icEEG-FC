#!/usr/bin/env python3
"""
Plot iEEG/fMRI connectome (nodes+edges) on FEAT highres anatomy.
*_bipolar_nodes_midpoints.csv is in example_func voxel coordinates

✅ 100% coordinate safety rule:
- The input bipolar midpoint CSV is in **example_func VOXEL space** (what you used for ROI extraction).
- To plot on **highres**, we convert:
      example_func (vox) -> example_func (mm) -> highres (vox)
  using NIfTI affines (no guessing, no sign hacks).

This guarantees the dots land correctly on highres, and you can still use the same FC matrices.

Written for Tara's ICE pipeline.
"""

import os
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =========================
# USER SETTINGS
# =========================
subj = "ICE070"

inputs_dir = f"/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/python_inputs/{subj}"
fmri_run = "Run1"
eeg_run  = "Run1"
# eeg_band = "delta"
# eeg_band = "theta"
# eeg_band = "alpha"
# eeg_band = "beta"
# eeg_band = "gammalow"
eeg_band = "gammahigh"


fmri_fc_csv = os.path.join(inputs_dir, f"{subj}_{fmri_run}_FMRI_FC_static_labeled.csv")
eeg_fc_csv  = os.path.join(inputs_dir, f"{subj}_{eeg_run}_EEG_FC_static_{eeg_band}_labeled.csv")

fmri_seed_csv = os.path.join(inputs_dir, f"{subj}_{fmri_run}_seedNames.csv")
eeg_seed_csv  = os.path.join(inputs_dir, f"{subj}_{eeg_run}_seedNames.csv")


# Bipolar midpoints you currently have (these are example_func VOXEL coords)
bip_csv  = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_bipolar_nodes_midpoints.csv"

# FEAT reg folder
reg_dir = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}/4_MRI/2_Functionals/2_PreProcessing/{fmri_run}.feat/reg"

# Backgrounds
func_nii   = os.path.join(reg_dir, "example_func.nii.gz")   # source space for coords
highres_nii= os.path.join(reg_dir, "highres.nii.gz")        # target space for plotting

out_dir = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_outputs_highresConnectome"
os.makedirs(out_dir, exist_ok=True)

# Plot knobs
NODE_SIZE = 40
EDGE_LW_MIN = 0.5
EDGE_LW_MAX = 5.0
TOP_PCT_EDGES = 15      # if set to 10 ==> show top 10% strongest |FC| (since threshold = 100-top_pct)
MAX_EDGES_CAP = 120
NEAR_Z_TOL = 20         # voxels around the slice
ONE_SLICE = True        # True => single axial slice at median z

# =========================
# HELPERS
# =========================
def norm(s: str) -> str:
    return str(s).strip().upper().replace("_", "-")

def load_fc(csv_path: str) -> pd.DataFrame:
    fc = pd.read_csv(csv_path, index_col=0)
    fc.index = fc.index.astype(str).map(norm)
    fc.columns = fc.columns.astype(str).map(norm)
    if fc.shape[0] == fc.shape[1]:
        fc = (fc + fc.T) / 2.0
    return fc

def load_bipolar_coords_vox(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "bipolar_name" not in df.columns:
        raise ValueError("bipolar midpoint CSV must have column 'bipolar_name'")
    for c in ["x", "y", "z"]:
        if c not in df.columns:
            raise ValueError(f"bipolar midpoint CSV must have columns x,y,z. Missing {c}")
    df["bipolar_name"] = df["bipolar_name"].map(norm)
    df = df.dropna(subset=["x","y","z"]).copy()
    return df

def vox_to_world(affine, ijk):
    """ijk (N,3) voxel -> xyz_mm (N,3)"""
    ijk = np.asarray(ijk, float)
    xyz1 = np.c_[ijk, np.ones((ijk.shape[0], 1))]
    xyz = (affine @ xyz1.T).T[:, :3]
    return xyz

def world_to_vox(affine, xyz_mm):
    """xyz_mm (N,3) -> ijk (N,3) voxel (float)"""
    invA = np.linalg.inv(affine)
    xyz_mm = np.asarray(xyz_mm, float)
    xyz1 = np.c_[xyz_mm, np.ones((xyz_mm.shape[0], 1))]
    ijk = (invA @ xyz1.T).T[:, :3]
    return ijk

def sanity_report_mapping(func_img, high_img, nodes_func_vox, nodes_high_vox, tag=""):
    fshp = func_img.shape[:3]
    hshp = high_img.shape[:3]

    inside_func = (
        (nodes_func_vox[:,0] >= 0) & (nodes_func_vox[:,0] < fshp[0]) &
        (nodes_func_vox[:,1] >= 0) & (nodes_func_vox[:,1] < fshp[1]) &
        (nodes_func_vox[:,2] >= 0) & (nodes_func_vox[:,2] < fshp[2])
    )
    inside_high = (
        (nodes_high_vox[:,0] >= 0) & (nodes_high_vox[:,0] < hshp[0]) &
        (nodes_high_vox[:,1] >= 0) & (nodes_high_vox[:,1] < hshp[1]) &
        (nodes_high_vox[:,2] >= 0) & (nodes_high_vox[:,2] < hshp[2])
    )

    print(f"\n[SANITY {tag}] func shape={fshp} axcodes={nib.aff2axcodes(func_img.affine)}")
    print(f"[SANITY {tag}] high shape={hshp} axcodes={nib.aff2axcodes(high_img.affine)}")
    print(f"[SANITY {tag}] inside func grid: {inside_func.sum()}/{len(nodes_func_vox)}")
    print(f"[SANITY {tag}] inside high grid: {inside_high.sum()}/{len(nodes_high_vox)}")

    print(f"[SANITY {tag}] func ijk ranges: "
          f"i[{nodes_func_vox[:,0].min():.1f},{nodes_func_vox[:,0].max():.1f}] "
          f"j[{nodes_func_vox[:,1].min():.1f},{nodes_func_vox[:,1].max():.1f}] "
          f"k[{nodes_func_vox[:,2].min():.1f},{nodes_func_vox[:,2].max():.1f}]")
    print(f"[SANITY {tag}] high ijk ranges: "
          f"i[{nodes_high_vox[:,0].min():.1f},{nodes_high_vox[:,0].max():.1f}] "
          f"j[{nodes_high_vox[:,1].min():.1f},{nodes_high_vox[:,1].max():.1f}] "
          f"k[{nodes_high_vox[:,2].min():.1f},{nodes_high_vox[:,2].max():.1f}]")

def build_edges(A: np.ndarray, top_pct=90, max_cap=120):
    """
    Keep strongest |FC| edges.
    top_pct=90 => threshold at 10th percentile? No:
      thr = percentile(scores, 100-top_pct) => percentile(scores, 10)
      keep scores >= 10th percentile => keeps top 90% (a lot).
    If you want "top 2%" edges, set top_pct=2.
    """
    n = A.shape[0]
    iu = np.triu_indices(n, k=1)
    vals = A[iu]
    scores = np.abs(vals)

    m = np.isfinite(scores) & (scores > 0)
    if not np.any(m):
        return []

    vals = vals[m]
    scores = scores[m]
    ii = iu[0][m]
    jj = iu[1][m]

    thr = np.percentile(scores, 100 - top_pct)
    keep = scores >= thr
    if not np.any(keep):
        return []

    ii, jj, vals, scores = ii[keep], jj[keep], vals[keep], scores[keep]
    order = np.argsort(scores)[::-1]
    ii, jj, vals, scores = ii[order], jj[order], vals[order], scores[order]

    if len(ii) > max_cap:
        ii, jj, vals, scores = ii[:max_cap], jj[:max_cap], vals[:max_cap], scores[:max_cap]

    return [(int(i), int(j), float(w), float(s)) for i, j, w, s in zip(ii, jj, vals, scores)]

def lw_from_scores(scores, lw_min=0.5, lw_max=5.0):
    scores = np.asarray(scores, float)
    if len(scores) == 0:
        return np.array([])
    smin, smax = scores.min(), scores.max()
    if smax == smin:
        return np.full_like(scores, (lw_min + lw_max) / 2.0)
    return lw_min + (scores - smin) * (lw_max - lw_min) / (smax - smin)

def pick_axial_slices(k_center, n_slices=5, step=6):
    half = n_slices // 2
    return [k_center + (i-half)*step for i in range(n_slices)]

def plot_connectome_on_highres(high_img, nodes_high_vox, node_names, edges, title, out_png):
    anat = high_img.get_fdata()
    shp = anat.shape[:3]

    ijk = np.rint(nodes_high_vox).astype(int)
    ijk[:, 0] = np.clip(ijk[:, 0], 0, shp[0]-1)
    ijk[:, 1] = np.clip(ijk[:, 1], 0, shp[1]-1)
    ijk[:, 2] = np.clip(ijk[:, 2], 0, shp[2]-1)

    k_med = int(np.median(ijk[:, 2]))
    if ONE_SLICE:
        ks = [k_med]
    else:
        ks = [k for k in pick_axial_slices(k_med, n_slices=5, step=6) if 0 <= k < shp[2]]

    scores = [s for (_, _, _, s) in edges]
    lws = lw_from_scores(scores, EDGE_LW_MIN, EDGE_LW_MAX)

    fig, axes = plt.subplots(1, len(ks), figsize=(4*len(ks), 5))
    if len(ks) == 1:
        axes = [axes]

    for ax, k in zip(axes, ks):
        ax.imshow(anat[:, :, k].T, origin="lower", cmap="gray")

        near = np.abs(ijk[:, 2] - k) <= NEAR_Z_TOL

        ax.scatter(
            ijk[near, 0], ijk[near, 1],
            s=NODE_SIZE,
            facecolors=(0.35, 0.70, 1.0, 0.70),
            edgecolors=(0.0, 0.35, 0.75),
            linewidths=1.2
        )

        for (idx, (i, j, w, s)) in enumerate(edges):
            if not (near[i] and near[j]):
                continue
            x1, y1 = ijk[i, 0], ijk[i, 1]
            x2, y2 = ijk[j, 0], ijk[j, 1]
            ls = "-" if w >= 0 else "--"
            ax.plot([x1, x2], [y1, y2], linewidth=lws[idx], linestyle=ls)

        ax.set_title(f"z={k}")
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close(fig)
    print("[OK] Saved:", out_png)

def load_seednames_csv(path: str):
    s = pd.read_csv(path, header=None).iloc[:,0].astype(str).map(norm).tolist()
    return s

fmri_seeds = load_seednames_csv(fmri_seed_csv)
eeg_seeds  = load_seednames_csv(eeg_seed_csv)

# sanity: FC labels should equal seed list (order matters if you expect it)
fc = load_fc(fmri_fc_csv)
if fc.index.tolist() != fmri_seeds:
    print("[WARN] fMRI seedNames order/contents != FC CSV labels.")
    # If you want to hard-fail instead:
    # raise ValueError("SeedNames CSV does not match FC CSV labels.")


# =========================
# MAIN
# =========================
for p in [func_nii, highres_nii, bip_csv, fmri_fc_csv, eeg_fc_csv, fmri_seed_csv, eeg_seed_csv]:
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing: {p}")


print("[INFO] Loading func:", func_nii)
func_img = nib.load(func_nii)
print("[INFO] Loading highres:", highres_nii)
high_img = nib.load(highres_nii)

coords_df = load_bipolar_coords_vox(bip_csv)
coords_df.to_csv(os.path.join(out_dir, f"{subj}_bipolar_midpoints_DEBUG.csv"), index=False)

coord_map_func_vox = {
    r["bipolar_name"]: (float(r["x"]), float(r["y"]), float(r["z"]))
    for _, r in coords_df.iterrows()
}

def run_one(modality_name, fc_csv, out_png):
    fc = load_fc(fc_csv)

    fc_nodes = list(fc.index)
    matched = [n for n in fc_nodes if n in coord_map_func_vox]
    missing = [n for n in fc_nodes if n not in coord_map_func_vox]

    print(f"\n[{modality_name}] FC nodes={len(fc_nodes)} matched coords={len(matched)} missing={len(missing)}")
    if missing:
        print("[MISSING examples]", missing[:30])

    if len(matched) < 2:
        raise RuntimeError(f"{modality_name}: too few matched nodes to plot.")

    # --- nodes in func voxel ---
    nodes_func_vox = np.array([coord_map_func_vox[n] for n in matched], float)

    # --- convert: func vox -> func mm -> highres vox ---
    nodes_mm = vox_to_world(func_img.affine, nodes_func_vox) # Convert those voxels to world (mm) using example_func affine
    nodes_high_vox = world_to_vox(high_img.affine, nodes_mm) # Convert those same mm points into highres voxel coordinates using highres affine

    sanity_report_mapping(func_img, high_img, nodes_func_vox, nodes_high_vox, tag=modality_name)

    # drop any nodes that fall outside highres grid (should be 0 drops if everything is correct)
    hshp = high_img.shape[:3]
    inside_high = (
        (nodes_high_vox[:,0] >= 0) & (nodes_high_vox[:,0] < hshp[0]) &
        (nodes_high_vox[:,1] >= 0) & (nodes_high_vox[:,1] < hshp[1]) &
        (nodes_high_vox[:,2] >= 0) & (nodes_high_vox[:,2] < hshp[2])
    )
    if not np.all(inside_high):
        bad = np.where(~inside_high)[0]
        print(f"[WARN] {modality_name}: dropping {len(bad)} nodes outside highres grid (check registration).")
        matched_keep = [n for i,n in enumerate(matched) if inside_high[i]]
        nodes_high_vox = nodes_high_vox[inside_high]
        matched = matched_keep

    # --- FC subset ---
    fc2 = fc.loc[matched, matched]
    A = fc2.to_numpy(float)
    np.fill_diagonal(A, 0.0)

    edges = build_edges(A, top_pct=TOP_PCT_EDGES, max_cap=MAX_EDGES_CAP)
    print(f"[{modality_name}] edges shown={len(edges)} (top {TOP_PCT_EDGES}% by |FC|, cap={MAX_EDGES_CAP})")

    plot_connectome_on_highres(
        high_img,
        nodes_high_vox,
        matched,
        edges,
        title=f"{subj} {modality_name} (HIGHRES axial) top {TOP_PCT_EDGES}% edges",
        out_png=out_png
    )

# fMRI
run_one(
    "fMRI static",
    fmri_fc_csv,
    os.path.join(out_dir, f"{subj}_fMRI_HIGHRES_axial.png")
)

# EEG
run_one(
    f"iEEG {eeg_band} static",
    eeg_fc_csv,
    os.path.join(out_dir, f"{subj}_iEEG_{eeg_band}_HIGHRES_axial.png")
)

print("\n[DONE]")
