#!/usr/bin/env python3

import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from ase import io
from copy import deepcopy

# --- CONFIGURATION ---
NEB_MEP_XYZ = 'neb_mep.xyz'      # Path to your final MEP xyz trajectory
NEB_ALL_XYZ = 'neb_all.xyz'      # Path to all NEB images (including prior iterations)
ERR_TOL = 0.01                  # RMS error tolerance in Angstrom
MAX_PCS = 6                     # Maximum allowed number of PCs
N_PRIOR = 4                      # Number of prior NEB iterations to validate on
PCA_OUT = 'network_parameters.npz'
INTERACTIVE = True  # Set to True for local/interactive use
"""
Notes:
1. All future workflow steps should be run on the aligned path created here
2. By default this accepts the lowest number of PCs that satisfies the rms training error specified;
   If you wish to use a specific number of PCs, set the tolerance to 0 and it will default to MAX_PCS
   or the max allowed number based on the number of images for training and the number of atoms in the system
   (roughly the total degrees of freedom).
Assumptions:
1. Each NEB iteration has the same number of images (and the same number as the final MEP)
2. All iterations are well-formed such that the total images in the ALL file is a multiple of that in the MEP
3. The atom indices are consistent in all images (and grouped by element, so they can be passed to QMCPACK later).
   Note that the indexing form you train on is the assumed form you will be running stalk on.
"""

# --- LOAD TRAJECTORIES ---
mep_images = io.read(NEB_MEP_XYZ + '@:')
all_images = io.read(NEB_ALL_XYZ + '@:')

n_atoms = mep_images[0].get_global_number_of_atoms()
traj_mep = np.array([img.positions.flatten() for img in mep_images])

# Determine how many images per iteration (assume all iterations have same count as MEP)
images_per_iter = len(traj_mep)
# Get the last N_PRIOR iterations (excluding the MEP itself)
traj_prior = []
if len(all_images) > images_per_iter:
    n_prior_total = min(N_PRIOR, (len(all_images) // images_per_iter) - 1)
    for i in range(n_prior_total):
        start = -(i+2)*images_per_iter
        end = -(i+1)*images_per_iter
        traj_prior.append(np.array([img.positions.flatten() for img in all_images[start:end]]))
    traj_prior = np.concatenate(traj_prior, axis=0)
else:
    traj_prior = None

# --- ALIGN (OPTIONAL, align all to MEP mean) ---
def kabsch_align(P, Q):
    C = np.dot(P.T, Q)
    V, S, W = np.linalg.svd(C)
    d = np.linalg.det(np.dot(V, W)) < 0.0
    if d:
        V[:, -1] *= -1
    U = np.dot(V, W)
    return U

def align_to_mean(traj, mean):
    n_images, n_flat = traj.shape
    n_atoms = n_flat // 3
    traj_reshaped = traj.reshape(n_images, n_atoms, 3)
    aligned = []
    for img in traj_reshaped:
        P = img - img.mean(axis=0)
        Q = mean - mean.mean(axis=0)
        U = kabsch_align(P, Q)
        aligned_img = np.dot(P, U) + mean.mean(axis=0)
        aligned.append(aligned_img.flatten())
    return np.array(aligned)

mep_mean = np.mean(traj_mep.reshape(len(traj_mep), n_atoms, 3), axis=0)
traj_mep_aligned = align_to_mean(traj_mep, mep_mean)
if traj_prior is not None:
    traj_prior_aligned = align_to_mean(traj_prior, mep_mean)

# --- SAVE ALIGNED XYZ ---
def get_aligned_xyz_name(xyz_path):
    base = os.path.basename(xyz_path)
    if '.' in base:
        base = '.'.join(base.split('.')[:-1])
    return base + '_aligned.xyz'

aligned_xyz_name = get_aligned_xyz_name(NEB_MEP_XYZ)
aligned_images = []
for orig_img, pos in zip(mep_images, traj_mep_aligned.reshape(-1, n_atoms, 3)):
    img = deepcopy(orig_img)
    img.positions = pos
    aligned_images.append(img)

io.write(aligned_xyz_name, aligned_images)
print(f"Aligned MEP written to {aligned_xyz_name}")

# --- STANDARDIZE (fit only on MEP) ---
scaler = StandardScaler(with_std=False)
traj_mep_scaled = scaler.fit_transform(traj_mep_aligned)
if traj_prior is not None:
    traj_prior_scaled = scaler.transform(traj_prior_aligned)

# --- PCA TRAINING (on MEP only) ---
max_possible_pcs = min(MAX_PCS, traj_mep_scaled.shape[0], traj_mep_scaled.shape[1])
pca = PCA(n_components=max_possible_pcs)
pca.fit(traj_mep_scaled)

# --- DIAGNOSTICS ---
def reconstruction_error(orig, recon, n_atoms):
    return np.sqrt(np.mean(np.sum((orig.reshape(-1, n_atoms, 3) - recon.reshape(-1, n_atoms, 3))**2, axis=2)))

print("\nPC | Train RMS error | Valid RMS error | Cumulative Var (%)")
print("---|-----------------|----------------|-------------------")
cum_var = 0
rms_train = []
rms_valid = []
for n in range(1, max_possible_pcs+1):
    pca_n = PCA(n_components=n)
    pca_n.fit(traj_mep_scaled)
    # Training error
    proj_train = pca_n.transform(traj_mep_scaled)
    recon_train = scaler.inverse_transform(pca_n.inverse_transform(proj_train))
    err_train = reconstruction_error(traj_mep_aligned, recon_train, n_atoms)
    rms_train.append(err_train)
    # Validation error
    if traj_prior is not None:
        proj_valid = pca_n.transform(traj_prior_scaled)
        recon_valid = scaler.inverse_transform(pca_n.inverse_transform(proj_valid))
        err_valid = reconstruction_error(traj_prior_aligned, recon_valid, n_atoms)
        rms_valid.append(err_valid)
    else:
        err_valid = float('nan')
    cum_var = np.sum(pca_n.explained_variance_ratio_)
    print(f"{n:2d} |     {err_train:8.6f}    |    {err_valid:8.6f}   |    {cum_var*100:8.6f}")

# --- SELECT NUMBER OF PCs ---
auto_n_pcs = next((i+1 for i, err in enumerate(rms_train) if err <= ERR_TOL), max_possible_pcs)
print(f"\nAuto-selected number of PCs: {auto_n_pcs} (Train RMS error = {rms_train[auto_n_pcs-1]:.4f} Å, threshold = {ERR_TOL})")

# --- INTERACTIVE SELECTION ---
if INTERACTIVE:
    try:
        user_input = input(f"Enter number of PCs to use (1-{max_possible_pcs}, Enter for auto): ")
        if user_input.strip():
            n_pcs = int(user_input)
            if not (1 <= n_pcs <= max_possible_pcs):
                print(f"Invalid input, using auto-selected value: {auto_n_pcs}")
                n_pcs = auto_n_pcs
        else:
            n_pcs = auto_n_pcs
    except Exception:
        n_pcs = auto_n_pcs
else:
    n_pcs = auto_n_pcs

print(f"Using {n_pcs} principal components.")

# --- FINAL PCA AND SAVE ---
pca_final = PCA(n_components=n_pcs)
pca_final.fit(traj_mep_scaled)
np.savez(PCA_OUT,
         components=pca_final.components_,
         explained_variance=pca_final.explained_variance_,
         explained_variance_ratio=pca_final.explained_variance_ratio_,
         pca_mean=pca_final.mean_ if hasattr(pca_final, 'mean_') else 0,
         scaler_mean=scaler.mean_,
         scaler_scale=scaler.scale_ if hasattr(scaler, 'scale_') else 1)

print(f"\nPCA model and scaler saved to {PCA_OUT}")