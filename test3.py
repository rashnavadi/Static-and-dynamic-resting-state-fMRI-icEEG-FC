import numpy as np
import pandas as pd
import nibabel as nib

subj = "ICE055"
anat_nii = f"/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/{subj}/4_MRI/1_Anatomicals/3D_Anat/{subj}_3Danat_brain.nii.gz"
bip_csv  = f"/Volumes/MIND/ICE/Tara/native_space_electrodes_coords/{subj}_bipolar_nodes_midpoints.csv"

img = nib.load(anat_nii)
A = img.affine
shp = img.shape[:3]
axcodes = nib.aff2axcodes(A)
print("[ANAT] axcodes:", axcodes, "shape:", shp)

# world-x at left-most voxel and right-most voxel (same y/z mid)
jmid, kmid = shp[1]//2, shp[2]//2
def vox2world(i,j,k):
    v = np.array([i,j,k,1.0])
    return (A @ v)[:3]

x0   = vox2world(0,        jmid, kmid)[0]
xend = vox2world(shp[0]-1, jmid, kmid)[0]
x_min, x_max = min(x0,xend), max(x0,xend)
x_mid = 0.5*(x_min + x_max)
print("[ANAT] world-x range:", x_min, "to", x_max, "midline x_mid≈", x_mid)

df = pd.read_csv(bip_csv)
df["bipolar_name"] = df["bipolar_name"].astype(str)

# classify by midline (NOT by sign!)
df["side_by_midline"] = np.where(df["x"] < x_mid, "side_A", "side_B")

print("\n[COORDS] x min/max:", df["x"].min(), df["x"].max())
print("[COORDS] counts by midline side:\n", df["side_by_midline"].value_counts())

print("\n[COORDS] first 30 with side:")
print(df[["bipolar_name","x","y","z","side_by_midline"]].head(30).to_string(index=False))
