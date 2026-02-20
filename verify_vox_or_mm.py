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

example_func_nii = "/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis/ICE013/4_MRI/2_Functionals/2_PreProcessing/Run1a.feat/reg/example_func.nii.gz"


img = nib.load(example_func_nii)
print(img.header.get_zooms())
print(img.affine)
