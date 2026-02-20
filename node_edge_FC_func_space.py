#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib

# =========================
# USER SETTINGS
# =========================
subj = "ICE055"

# anat_nii = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}/4_MRI/1_Anatomicals/3D_Anat/{subj}_3Danat_brain.nii.gz"

bip_csv  = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_bipolar_nodes_midpoints.csv"

inputs_dir = f"/Volumes/MIND/ICE/Tara/Kristina_paper/FC_rebuild/RESULTS_REPRESENTATIVE/python_inputs/{subj}"
fmri_run = "Run1a"
eeg_run  = "Run1a"
eeg_band = "delta"
func_nii = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}/4_MRI/2_Functionals/2_PreProcessing/{fmri_run}.feat/reg/example_func.nii.gz"


fmri_fc_csv = os.path.join(inputs_dir, f"{subj}_{fmri_run}_FMRI_FC_static_labeled.csv")
eeg_fc_csv  = os.path.join(inputs_dir, f"{subj}_{eeg_run}_EEG_FC_static_{eeg_band}_labeled.csv")

out_dir = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_outputs_nativeConnectome"
os.makedirs(out_dir, exist_ok=True)

def sanity_check_coords_in_anat(func_nii, xyz_mm, tag="coords"):
    img = nib.load(func_nii)
    A = img.affine
    shp = img.shape[:3]

    # voxel coords
    invA = np.linalg.inv(A)
    ijk = (invA @ np.c_[xyz_mm, np.ones(len(xyz_mm))].T).T[:, :3]

    inside = (
        (ijk[:,0] >= 0) & (ijk[:,0] < shp[0]) &
        (ijk[:,1] >= 0) & (ijk[:,1] < shp[1]) &
        (ijk[:,2] >= 0) & (ijk[:,2] < shp[2])
    )

    axcodes = nib.aff2axcodes(A)
    print(f"[SANITY] {tag}: {inside.sum()}/{len(xyz_mm)} points fall INSIDE anat voxel grid")
    print(f"[SANITY] anat orientation: {axcodes}, shape={shp}")
    print(f"[SANITY] ijk ranges: i[{ijk[:,0].min():.1f},{ijk[:,0].max():.1f}] "
          f"j[{ijk[:,1].min():.1f},{ijk[:,1].max():.1f}] "
          f"k[{ijk[:,2].min():.1f},{ijk[:,2].max():.1f}]")

    return inside

# converts mm → voxel
def world_to_voxel(affine, xyz_mm):
    invA = np.linalg.inv(affine)
    xyz1 = np.c_[xyz_mm, np.ones((xyz_mm.shape[0], 1))]
    ijk  = (invA @ xyz1.T).T[:, :3]
    return ijk

def report_lr(anat_img, nodes_vox, node_names):
    axcodes = nib.aff2axcodes(anat_img.affine)  # e.g. ('L','A','S') or ('R','A','S')
    shp = anat_img.shape[:3]
    i_mid = (shp[0] - 1) / 2.0

    ijk = world_to_voxel(anat_img.affine, nodes_vox)
    i = ijk[:, 0]

    # side by voxel i relative to mid
    left_vox  = i > i_mid
    right_vox = i <= i_mid

    # interpret which side is "Left" depending on axis code
    # If first axis is 'L', increasing i goes Left. If 'R', increasing i goes Right.
    if axcodes[0] == 'L':
        left_mask, right_mask = left_vox, right_vox
    elif axcodes[0] == 'R':
        left_mask, right_mask = right_vox, left_vox
    else:
        print("[WARN] Unexpected x-axis code:", axcodes[0], "Cannot label Left/Right reliably.")
        left_mask, right_mask = left_vox, right_vox

    print("\n[ORIENT]", axcodes, "shape=", shp, "i_mid=", i_mid)
    print("[HEMISPHERE COUNTS] Left =", int(left_mask.sum()), "Right =", int(right_mask.sum()))

    # show a few examples
    print("\nExamples classified LEFT:")
    for n in np.array(node_names)[left_mask][:8]:
        print(" ", n)
    print("Examples classified RIGHT:")
    for n in np.array(node_names)[right_mask][:8]:
        print(" ", n)

# usage (after you build nodes_vox + node_names):
# report_lr(anat_img, nodes_vox, matched_names)



# Plot knobs
NODE_SIZE = 40          # scatter marker size
EDGE_LW_MIN = 0.5       # min line width
EDGE_LW_MAX = 5.0       # max line width
TOP_PCT_EDGES = 90       # show top 2% edges by |FC|
MAX_EDGES_CAP = 120     # hard cap for readability

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

def load_bipolar_coords(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "bipolar_name" not in df.columns:
        raise ValueError("bipolar midpoint CSV must have column 'bipolar_name'")
    for c in ["x", "y", "z"]:
        if c not in df.columns:
            raise ValueError(f"bipolar midpoint CSV must have columns x,y,z. Missing {c}")
    df["bipolar_name"] = df["bipolar_name"].map(norm)
    df = df.dropna(subset=["x","y","z"]).copy()
    return df

def world_to_voxel(affine, xyz_mm):
    invA = np.linalg.inv(affine)
    xyz1 = np.c_[xyz_mm, np.ones((xyz_mm.shape[0], 1))]
    ijk  = (invA @ xyz1.T).T[:, :3]
    return ijk

def build_edges(A: np.ndarray, top_pct=2, max_cap=120):
    """
    Returns list of (i, j, w, s) where:
      w = signed FC
      s = abs(FC) score used for ranking
    """
    n = A.shape[0]
    iu = np.triu_indices(n, k=1)
    vals = A[iu]
    scores = np.abs(vals)

    # drop zeros / nans
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

    edges = [(int(i), int(j), float(w), float(s)) for i, j, w, s in zip(ii, jj, vals, scores)]
    return edges

def lw_from_scores(scores, lw_min=0.5, lw_max=5.0):
    scores = np.asarray(scores, float)
    if len(scores) == 0:
        return np.array([])
    smin, smax = scores.min(), scores.max()
    if smax == smin:
        return np.full_like(scores, (lw_min + lw_max)/2)
    return lw_min + (scores - smin) * (lw_max - lw_min) / (smax - smin)

def pick_axial_slices(k_center, n_slices=5, step=6):
    half = n_slices // 2
    return [k_center + (i-half)*step for i in range(n_slices)]

def plot_native_connectome(anat_img, nodes_vox, node_names, edges, title, out_png):
    anat = anat_img.get_fdata()
    A = anat_img.affine
    shp = anat.shape[:3]

    # node voxel coords (for plotting on slices)
    # ijk = world_to_voxel(A, nodes_vox)
    # ijk_round = np.rint(ijk).astype(int)
    ijk_round = np.rint(nodes_vox).astype(int)   # already voxel coords
    ijk_round[:, 0] = np.clip(ijk_round[:, 0], 0, shp[0]-1)
    ijk_round[:, 1] = np.clip(ijk_round[:, 1], 0, shp[1]-1)
    ijk_round[:, 2] = np.clip(ijk_round[:, 2], 0, shp[2]-1)

    # pick axial slices around median z
    k_med = int(np.median(ijk_round[:, 2]))
    ks = [k for k in pick_axial_slices(k_med, n_slices=5, step=6) if 0 <= k < shp[2]]

    # line widths
    scores = [s for (_, _, _, s) in edges]
    lws = lw_from_scores(scores, EDGE_LW_MIN, EDGE_LW_MAX)

    fig, axes = plt.subplots(1, len(ks), figsize=(4*len(ks), 5))
    if len(ks) == 1:
        axes = [axes]

    for ax, k in zip(axes, ks):
        # axial slice: (x,y,k) -> show as y-by-x with transpose
        ax.imshow(anat[:, :, k].T, origin="lower", cmap="gray")
        
        # ax.invert_xaxis() # match radiological display


        # nodes near this slice (within +-2 vox)
        near = np.abs(ijk_round[:, 2] - k) <= 20
        # near = np.ones(len(ijk_round), dtype=bool)   # show all nodes on all slices

        ax.scatter(ijk_round[near, 0], ijk_round[near, 1], s=NODE_SIZE)

        # edges: draw if both endpoints are near the slice
        for (idx, (i, j, w, s)) in enumerate(edges):
            if not (near[i] and near[j]):
                continue
            x1, y1 = ijk_round[i, 0], ijk_round[i, 1]
            x2, y2 = ijk_round[j, 0], ijk_round[j, 1]
            # sign: positive solid, negative dashed (simple visual cue)
            ls = "-" if w >= 0 else "--"
            ax.plot([x1, x2], [y1, y2], linewidth=lws[idx], linestyle=ls)

        ax.set_title(f"z={k}")
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close(fig)
    print("[OK] Saved:", out_png)

# =========================
# MAIN
# =========================
print("Loading anatomy:", func_nii)
anat_img = nib.load(func_nii)

coords_df = load_bipolar_coords(bip_csv)
coords_df.to_csv(os.path.join(out_dir, f"{subj}_bipolar_midpoints_DEBUG.csv"), index=False)

coord_map = {r["bipolar_name"]: (float(r["x"]), float(r["y"]), float(r["z"]))
             for _, r in coords_df.iterrows()}


print("anat shape:", anat_img.shape)
print("coord x range:", coords_df["x"].min(), coords_df["x"].max())

def run_one(modality_name, fc_csv, out_png):
    fc = load_fc(fc_csv)

    fc_nodes = list(fc.index)
    matched = [n for n in fc_nodes if n in coord_map]
    missing = [n for n in fc_nodes if n not in coord_map]

    print(f"\n[{modality_name}] FC nodes={len(fc_nodes)} matched coords={len(matched)} missing={len(missing)}")
    if missing:
        print("[MISSING examples]", missing[:30])

    print("[MATCHED nodes + coords]")
    for n in matched:
        x, y, z = coord_map[n]
        print(f"  {n:15s}  x={x:9.3f} y={y:9.3f} z={z:9.3f}")

    if len(matched) < 2:
        raise RuntimeError(f"{modality_name}: too few matched nodes to plot.")

    fc2 = fc.loc[matched, matched]
    A = fc2.to_numpy(float)
    np.fill_diagonal(A, 0.0)

    edges = build_edges(A, top_pct=TOP_PCT_EDGES, max_cap=MAX_EDGES_CAP)
    print(f"[{modality_name}] edges shown={len(edges)} (top {TOP_PCT_EDGES}% by |FC|, cap={MAX_EDGES_CAP})")

    nodes_vox = np.array([coord_map[n] for n in matched], float)
    plot_native_connectome(
        anat_img,
        nodes_vox,
        matched,
        edges,
        title=f"{subj} {modality_name} (native axial) top {TOP_PCT_EDGES}% edges",
        out_png=out_png
    )

# fMRI
run_one(
    "fMRI static",
    fmri_fc_csv,
    os.path.join(out_dir, f"{subj}_fMRI_native_axial.png")
)

# EEG
run_one(
    f"iEEG {eeg_band} static",
    eeg_fc_csv,
    os.path.join(out_dir, f"{subj}_iEEG_{eeg_band}_native_axial.png")
)

print("\n[DONE]")
