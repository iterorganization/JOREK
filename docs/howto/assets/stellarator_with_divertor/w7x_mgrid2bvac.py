"""
Append vacuum field from mgrid_w7x.nc to gvec2jorek.dat (use_ext_field must be .true.).

Pipeline:
  1. Read existing gvec2jorek.dat -> grid geometry R(s,theta,phi), Z(s,theta,phi)
  2. Read mgrid_w7x.nc -> B_vac(R,Z,phi) per coil group on cylindrical grid
  3. Sum coil contributions with EXTCUR weights
  4. For each (s_i, theta_j), evaluate R,Z at N_phi_sample toroidal angles
  5. Interpolate mgrid B_vac at those physical (R,Z,phi) positions
  6. Toroidal Fourier transform -> mode coefficients (n_modes)
  7. Finite-difference derivatives: _s, _t, _st
  8. Append 12 B_vac blocks to output gvec2jorek.dat

Usage:
  python mgrid2bvac.py

Requirements: numpy, scipy, netCDF4
"""

import sys
import shutil
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import netCDF4 as nc
from helper_gvec2jorek import read_gvec_header, read_gvec_blocks, write_gvec_block

# =============================================================================
# Configuration

# TO DO: add actual link to files when become available

# Input files - specify path where located
GVEC2JOREK_IN = Path("gvec2jorek.dat")
MGRID_FILE    = Path("mgrid_w7x.nc")
GVEC2JOREK_OUT = Path("gvec2jorek_with_bvac.dat") # copy of input + B_vac blocks appended

# Coil currents [A] for each coil group in mgrid (7 groups for W7-X)
# These multiply the B-per-unit-current stored in mgrid.
# W7-X EIM standard config -- CALIBRATED to match F0/R = -2.42 T at R=6.0
# Uniform estimate for modular coils, trim coils off.
# One has to verify with actual VMEC input extcur values for the specific equilibrium!
EXTCUR = [13470.0, 13470.0, 13470.0, 13470.0, 13470.0, 0.0, 0.0]

# Number of toroidal sample points for Fourier transform (per field period)
N_PHI_SAMPLE = 36 

# =============================================================================
# I/O FUNCTIONS
# =============================================================================

def evaluate_fourier(mode_data, phi_arr, nfp):
    """Evaluate toroidal Fourier modes at given phi angles.

    Args:
        mode_data: (n_theta, n_rad, n_modes) Fourier coefficients
        phi_arr: (N_phi,) toroidal angles [rad]
        nfp: number of field periods

    Returns:
        (n_theta, n_rad, N_phi) real-space values
    """
    n_theta, n_rad, n_modes = mode_data.shape
    n_max = (n_modes - 1) // 2
    N_phi = len(phi_arr)

    sin_n = np.arange(1, n_max + 1)
    cos_n = np.arange(0, n_max + 1)

    # (N_phi, n_sin) and (N_phi, n_cos)
    sin_vals = np.sin(nfp * sin_n[None, :] * phi_arr[:, None])
    cos_vals = np.cos(nfp * cos_n[None, :] * phi_arr[:, None])

    sin_coeffs = mode_data[:, :, :n_max]        # (n_theta, n_rad, n_max)
    cos_coeffs = mode_data[:, :, n_max:n_modes]  # (n_theta, n_rad, n_max+1)

    # GVEC convention: negative sine (right-hand coord system)
    result = (-np.einsum('ijk,lk->ijl', sin_coeffs, sin_vals) +
               np.einsum('ijk,lk->ijl', cos_coeffs, cos_vals))
    return result


def fourier_transform(real_data, phi_arr, nfp, n_max):
    """Toroidal Fourier transform: real space -> mode coefficients.

    Args:
        real_data: (n_theta, n_rad, N_phi) real-space values
        phi_arr: (N_phi,) toroidal angles [rad]
        nfp: number of field periods
        n_max: maximum toroidal mode number

    Returns:
        (n_theta, n_rad, 2*n_max+1) mode coefficients
        Layout: [sin(1)..sin(n_max), cos(0)..cos(n_max)]
    """
    n_theta, n_rad, N_phi = real_data.shape
    n_modes = 2 * n_max + 1

    sin_n = np.arange(1, n_max + 1)
    cos_n = np.arange(0, n_max + 1)

    sin_vals = np.sin(nfp * sin_n[None, :] * phi_arr[:, None])  # (N_phi, n_max)
    cos_vals = np.cos(nfp * cos_n[None, :] * phi_arr[:, None])  # (N_phi, n_max+1)

    # Least-squares fit: build basis matrix (N_phi, n_modes)
    # Columns: -sin(n*nfp*phi) for n=1..n_max, cos(n*nfp*phi) for n=0..n_max
    basis = np.zeros((N_phi, n_modes))
    basis[:, :n_max] = -sin_vals      # negative sin convention
    basis[:, n_max:] = cos_vals

    # Solve for each (theta, rad) point
    # real_data reshaped: (n_theta*n_rad, N_phi)
    data_flat = real_data.reshape(-1, N_phi)  # (n_theta*n_rad, N_phi)

    # Least squares: basis^T basis * coeffs = basis^T data
    # Use np.linalg.lstsq for robustness
    coeffs_flat, _, _, _ = np.linalg.lstsq(basis, data_flat.T, rcond=None)
    # coeffs_flat: (n_modes, n_theta*n_rad)

    coeffs = coeffs_flat.T.reshape(n_theta, n_rad, n_modes)
    return coeffs

# =============================================================================
# MGRID READING
# =============================================================================

def read_mgrid(filepath, extcur):
    """Read mgrid_w7x.nc and sum coil contributions.

    Returns:
        R_grid: (nR,) array of R values [m]
        Z_grid: (nZ,) array of Z values [m]
        phi_grid: (nphi,) array of phi values [rad] (one field period)
        B_R, B_Z, B_phi: (nphi, nZ, nR) total vacuum field [T]
    """
    ds = nc.Dataset(str(filepath), 'r')

    nR = int(ds.variables['ir'][...])
    nZ = int(ds.variables['jz'][...])
    nphi = int(ds.variables['kp'][...])
    nfp = int(ds.variables['nfp'][...])
    nextcur = int(ds.variables['nextcur'][...])
    Rmin = float(ds.variables['rmin'][...])
    Rmax = float(ds.variables['rmax'][...])
    Zmin = float(ds.variables['zmin'][...])
    Zmax = float(ds.variables['zmax'][...])
    raw_cur = ds.variables['raw_coil_cur'][:]

    R_grid = np.linspace(Rmin, Rmax, nR)
    Z_grid = np.linspace(Zmin, Zmax, nZ)
    phi_grid = np.linspace(0, 2*np.pi/nfp, nphi, endpoint=False)

    print(f"mgrid: {nR}x{nZ}x{nphi} grid, {nextcur} coil groups, nfp={nfp}")
    print(f"  R=[{Rmin:.2f}, {Rmax:.2f}] m, Z=[{Zmin:.2f}, {Zmax:.2f}] m")

    B_R = np.zeros((nphi, nZ, nR))
    B_Z = np.zeros((nphi, nZ, nR))
    B_phi = np.zeros((nphi, nZ, nR))

    for ic in range(nextcur):
        cur = extcur[ic]
        br = ds.variables[f'br_{ic+1:03d}'][:]  # (nphi, nZ, nR)
        bz = ds.variables[f'bz_{ic+1:03d}'][:]
        bp = ds.variables[f'bp_{ic+1:03d}'][:]
        B_R += cur * br
        B_Z += cur * bz
        B_phi += cur * bp
        print(f"  Coil {ic+1}: extcur={cur:.1f} A, raw_cur={raw_cur[ic]:.0f}")

    ds.close()

    B_mag = np.sqrt(B_R**2 + B_Z**2 + B_phi**2)
    midR = nR // 2
    midZ = nZ // 2
    print(f"  |B| at grid center (phi=0): {B_mag[0, midZ, midR]:.4f} T")
    print(f"  |B| range: [{B_mag.min():.4f}, {B_mag.max():.4f}] T")

    return R_grid, Z_grid, phi_grid, B_R, B_Z, B_phi, nfp


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    print("\nmgrid2bvac: Append vacuum field to gvec2jorek.dat")
    print("=" * 70)

    # --- Step 1: Read gvec2jorek.dat ---
    print("\n[1] Reading gvec2jorek.dat...")
    hdr = read_gvec_header(GVEC2JOREK_IN)
    n_rad = hdr['n_rad']
    n_theta = hdr['n_theta']
    n_modes = hdr['n_modes']
    n_max = hdr['n_max']
    nfp = hdr['nfp']
    print(f"  Grid: n_rad={n_rad}, n_theta={n_theta}, n_modes={n_modes}, n_max={n_max}, nfp={nfp}")

    blocks, block_names = read_gvec_blocks(GVEC2JOREK_IN, n_theta, n_rad, n_modes)
    print(f"  Read {len(blocks)} data blocks")

    # Check required geometry blocks exist
    for name in ['R', 'R_s', 'R_t', 'R_st', 'Z', 'Z_s', 'Z_t', 'Z_st']:
        if name not in blocks:
            print(f"  ERROR: Missing block '{name}' in gvec2jorek.dat")
            sys.exit(1)
    print("  Geometry blocks (R, Z) present")

    # --- Step 2: Read mgrid ---
    print("\n[2] Reading mgrid_w7x.nc...")
    R_mg, Z_mg, phi_mg, BR_mg, BZ_mg, Bphi_mg, nfp_mg = read_mgrid(MGRID_FILE, EXTCUR)
    if nfp_mg != nfp:
        print(f"  ERROR: nfp mismatch: mgrid={nfp_mg}, gvec={nfp}")
        sys.exit(1)

    # Build interpolators (phi, Z, R) -> B
    print("\n[3] Building RegularGridInterpolators...")
    interp_BR   = RegularGridInterpolator((phi_mg, Z_mg, R_mg), BR_mg,
                                          method='linear', bounds_error=False, fill_value=None)
    interp_BZ   = RegularGridInterpolator((phi_mg, Z_mg, R_mg), BZ_mg,
                                          method='linear', bounds_error=False, fill_value=None)
    interp_Bphi = RegularGridInterpolator((phi_mg, Z_mg, R_mg), Bphi_mg,
                                          method='linear', bounds_error=False, fill_value=None)
    print("  Done (linear interpolation, extrapolate at edges)")

    # --- Step 3: Evaluate geometry at sample toroidal angles ---
    print(f"\n[4] Evaluating grid geometry at {N_PHI_SAMPLE} toroidal planes...")
    phi_sample = np.linspace(0, 2*np.pi/nfp, N_PHI_SAMPLE, endpoint=False)

    # R(s, theta, phi) and Z(s, theta, phi): evaluate from Fourier modes
    R_real = evaluate_fourier(blocks['R'], phi_sample, nfp)     # (n_theta, n_rad, N_phi)
    Z_real = evaluate_fourier(blocks['Z'], phi_sample, nfp)     # (n_theta, n_rad, N_phi)

    print(f"  R range: [{R_real.min():.3f}, {R_real.max():.3f}] m")
    print(f"  Z range: [{Z_real.min():.3f}, {Z_real.max():.3f}] m")

    # Check if all points fall within mgrid domain
    R_min_mg, R_max_mg = R_mg[0], R_mg[-1]
    Z_min_mg, Z_max_mg = Z_mg[0], Z_mg[-1]
    n_outside_R = np.sum((R_real < R_min_mg) | (R_real > R_max_mg))
    n_outside_Z = np.sum((Z_real < Z_min_mg) | (Z_real > Z_max_mg))
    n_total = R_real.size
    if n_outside_R > 0 or n_outside_Z > 0:
        print(f"  WARNING: {n_outside_R} R-points outside [{R_min_mg}, {R_max_mg}]")
        print(f"  WARNING: {n_outside_Z} Z-points outside [{Z_min_mg}, {Z_max_mg}]")
        print(f"  ({n_outside_R + n_outside_Z}/{n_total} total outside, will extrapolate)")
    else:
        print("  All points within mgrid domain")

    # --- Step 4: Interpolate vacuum field ---
    print("\n[5] Interpolating vacuum field at JOREK grid points...")

    # Remember coordinate convention mapping.
    # gvec2jorek.dat Fourier modes use GVEC/JOREK left-handed toroidal angle
    # (zeta), while mgrid uses standard RH cylindrical phi.
    # Relation: phi_cylindrical = -zeta_GVEC = (2*pi/nfp - zeta) mod (2*pi/nfp)
    # The geometry evaluation above gives (R, Z) at zeta = phi_sample.
    # These physical points sit at cylindrical phi = 2*pi/nfp - phi_sample.
    # (Same convention as FLIP_DIVERTOR_PHI in w7x analysis scripts.)
    phi_cyl = (2 * np.pi / nfp - phi_sample) % (2 * np.pi / nfp)
    print(f"  Phi convention: zeta_GVEC -> phi_cyl = 2pi/nfp - zeta")

    # Build interpolation points: (N_total, 3) with (phi_cyl, Z, R)
    points = np.stack([
        np.broadcast_to(phi_cyl[None, None, :], R_real.shape).ravel(),
        Z_real.ravel(),
        R_real.ravel()
    ], axis=-1)

    Bvac_R_real   = interp_BR(points).reshape(n_theta, n_rad, N_PHI_SAMPLE)
    Bvac_Z_real   = interp_BZ(points).reshape(n_theta, n_rad, N_PHI_SAMPLE)
    Bvac_phi_real = interp_Bphi(points).reshape(n_theta, n_rad, N_PHI_SAMPLE)

    # Remember sign convention for toroidal component.
    # mgrid uses standard RH cylindrical (phi counter-clockwise from above).
    # GVEC/JOREK use LH system: (x,y,z)=(Rcos(zeta),-Rsin(zeta),Z), so zeta = -phi.
    # B_R and B_Z are the same in both conventions.
    # B_phi_GVEC = -B_phi_cylindrical (toroidal direction reversed).
    Bvac_phi_real = -Bvac_phi_real

    Bvac_mag = np.sqrt(Bvac_R_real**2 + Bvac_Z_real**2 + Bvac_phi_real**2)
    print(f"  |B_vac| range: [{Bvac_mag.min():.4f}, {Bvac_mag.max():.4f}] T")
    print(f"  |B_vac| at axis (s=0, theta=0, phi=0): {Bvac_mag[0,0,0]:.4f} T")

    # --- Step 5: Toroidal Fourier transform ---
    print(f"\n[6] Fourier transforming (n_max={n_max}, n_modes={n_modes})...")
    Bvac_R_four   = fourier_transform(Bvac_R_real,   phi_sample, nfp, n_max)
    Bvac_Z_four   = fourier_transform(Bvac_Z_real,   phi_sample, nfp, n_max)
    Bvac_phi_four = fourier_transform(Bvac_phi_real, phi_sample, nfp, n_max)
    print(f"  Mode coefficient arrays: {Bvac_R_four.shape}")

    # Verify roundtrip
    Bvac_R_check = evaluate_fourier(Bvac_R_four, phi_sample, nfp)
    err = np.max(np.abs(Bvac_R_check - Bvac_R_real))
    print(f"  Fourier roundtrip max error (B_R): {err:.2e}")

    # Per-radial-layer errors (important: large errors near s=1 from coil proximity)
    err_per_s = np.max(np.abs(Bvac_R_check - Bvac_R_real), axis=(0,2))
    print(f"  Roundtrip error at s=0.0: {err_per_s[0]:.2e} T")
    print(f"  Roundtrip error at s=0.5: {err_per_s[n_rad//2]:.2e} T")
    print(f"  Roundtrip error at s=1.0: {err_per_s[-1]:.2e} T")
    print(f"  (Large errors near s=1 expected from coil proximity -- high-n content filtered)")

    # Compare with GVEC equilibrium B field as sanity check
    B_phi_gvec = evaluate_fourier(blocks['B_phi'], phi_sample, nfp)
    # GVEC stores total B (vacuum + plasma), ours is vacuum only
    # At s=0, plasma current contribution is small -> should be close
    diff_axis = np.abs(Bvac_phi_real[0, 0, 0] - B_phi_gvec[0, 0, 0])
    print(f"  B_vac_phi vs GVEC B_phi at axis: {Bvac_phi_real[0,0,0]:.4f} vs {B_phi_gvec[0,0,0]:.4f} (diff={diff_axis:.4f})")

    # --- Step 6: Compute derivatives via finite differences ---
    print("\n[7] Computing s and theta derivatives (2nd order central FD)...")
    ds = 1.0 / (n_rad - 1)   # uniform spacing in s
    dt = 2 * np.pi / n_theta  # uniform spacing in theta

    def fd_s(arr):
        """2nd-order central finite difference in s (axis=1), with
        forward/backward at boundaries."""
        out = np.zeros_like(arr)
        # Interior: central
        out[:, 1:-1, :] = (arr[:, 2:, :] - arr[:, :-2, :]) / (2 * ds)
        # s=0 boundary: forward
        out[:, 0, :] = (-3*arr[:, 0, :] + 4*arr[:, 1, :] - arr[:, 2, :]) / (2 * ds)
        # s=1 boundary: backward
        out[:, -1, :] = (3*arr[:, -1, :] - 4*arr[:, -2, :] + arr[:, -3, :]) / (2 * ds)
        return out

    def fd_t(arr):
        """2nd-order central finite difference in theta (axis=0), periodic."""
        out = np.zeros_like(arr)
        out[1:-1, :, :] = (arr[2:, :, :] - arr[:-2, :, :]) / (2 * dt)
        # Periodic boundary
        out[0, :, :]  = (arr[1, :, :] - arr[-1, :, :]) / (2 * dt)
        out[-1, :, :] = (arr[0, :, :] - arr[-2, :, :]) / (2 * dt)
        return out

    # Derivatives
    Bvac_R_s    = fd_s(Bvac_R_four)
    Bvac_R_t    = fd_t(Bvac_R_four)
    Bvac_R_st_s_t = fd_s(Bvac_R_t)
    Bvac_R_st_t_s = fd_t(Bvac_R_s)
    Bvac_R_st   = 0.5 * (Bvac_R_st_s_t + Bvac_R_st_t_s)

    Bvac_Z_s    = fd_s(Bvac_Z_four)
    Bvac_Z_t    = fd_t(Bvac_Z_four)
    Bvac_Z_st_s_t = fd_s(Bvac_Z_t)
    Bvac_Z_st_t_s = fd_t(Bvac_Z_s)
    Bvac_Z_st   = 0.5 * (Bvac_Z_st_s_t + Bvac_Z_st_t_s)

    Bvac_phi_s  = fd_s(Bvac_phi_four)
    Bvac_phi_t  = fd_t(Bvac_phi_four)
    Bvac_phi_st_s_t = fd_s(Bvac_phi_t)
    Bvac_phi_st_t_s = fd_t(Bvac_phi_s)
    Bvac_phi_st = 0.5 * (Bvac_phi_st_s_t + Bvac_phi_st_t_s)

    # Report discrete mixed-derivative commutator before symmetrization.
    comm_R = np.max(np.abs(Bvac_R_st_s_t - Bvac_R_st_t_s))
    comm_Z = np.max(np.abs(Bvac_Z_st_s_t - Bvac_Z_st_t_s))
    comm_P = np.max(np.abs(Bvac_phi_st_s_t - Bvac_phi_st_t_s))
    print(f"  Mixed-derivative commutator max |d_s d_t - d_t d_s|:")
    print(f"    R:   {comm_R:.3e}")
    print(f"    Z:   {comm_Z:.3e}")
    print(f"    phi: {comm_P:.3e}")

    print("  Done")

    # --- Step 7: Write augmented gvec2jorek.dat ---
    print(f"\n[8] Writing output: {GVEC2JOREK_OUT}")

    # Copy original file
    shutil.copy2(GVEC2JOREK_IN, GVEC2JOREK_OUT)

    # Append B_vac blocks
    with open(GVEC2JOREK_OUT, 'a') as f:
        bvac_blocks = [
            ('B_vac_R',       Bvac_R_four),
            ('B_vac_R_s',     Bvac_R_s),
            ('B_vac_R_t',     Bvac_R_t),
            ('B_vac_R_st',    Bvac_R_st),
            ('B_vac_Z',       Bvac_Z_four),
            ('B_vac_Z_s',     Bvac_Z_s),
            ('B_vac_Z_t',     Bvac_Z_t),
            ('B_vac_Z_st',    Bvac_Z_st),
            ('B_vac_phi',     Bvac_phi_four),
            ('B_vac_phi_s',   Bvac_phi_s),
            ('B_vac_phi_t',   Bvac_phi_t),
            ('B_vac_phi_st',  Bvac_phi_st),
        ]
        for name, data in bvac_blocks:
            write_gvec_block(f, data, name)
            print(f"  Wrote block: {name}")

if __name__ == '__main__':
    main()
