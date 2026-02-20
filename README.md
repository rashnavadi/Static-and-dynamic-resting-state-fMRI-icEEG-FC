Static and Dynamic Resting-State fMRI–icEEG Functional Connectivity

This repository contains the full analysis pipeline for investigating static and dynamic cross-modal coupling between resting-state fMRI functional connectivity (FC) and intracranial EEG (icEEG) functional connectivity.

The workflow progresses from signal extraction to hierarchical linear mixed-effects modeling.

Pipeline Overview (Steps 00–11)

Steps 00–11 implement the complete analysis pipeline:

fMRI time series extraction

Static fMRI FC computation

icEEG Hilbert-based band-limited power (BLP) FC computation

Window segmentation for dynamic FC estimation

Master table construction (subject–run–edge–window level)

True static pair-level FC derivation

Static linear mixed-effects modeling (pair-level)

Dynamic linear mixed-effects modeling (window-level)

Spatial distance modeling

Motion covariate integration

Reliability analysis (intraclass correlation; ICC)

The pipeline enables:

Dissociation of broadband static coupling from frequency-specific dynamic deviations

Population-level inference using hierarchical modeling

Edge-level modeling across subjects, runs, and windows

Statistical Framework

All statistical analyses are performed using hierarchical linear mixed-effects models (LMMs) to account for:

Subject-level variability

Run-level structure

Edge-level dependence

Both static (pair-level) and dynamic (window-level) models are implemented.

Repository Structure
Step00_...
Step01_...
...
Step11_...

Each step is modular and can be run independently after required inputs are generated.

Notes

Designed for simultaneous icEEG–fMRI datasets

Supports band-specific electrophysiological analyses

Includes motion and spatial distance covariates
