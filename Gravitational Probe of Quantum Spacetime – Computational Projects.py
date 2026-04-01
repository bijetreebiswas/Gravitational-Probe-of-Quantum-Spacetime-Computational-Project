# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 09:59:52 2026

@author: bijet
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d

# ============================================================
# 1.  POTENTIAL DEFINITIONS (from the paper)
# ============================================================
R = 2.0               # Schwarzschild radius (M=1)

def V_RW(r, l):
    """Classical Regge-Wheeler potential."""
    return (r - R) * (l*(l+1)*r - 3*R) / r**4

def V_q(r, qm, l):
    """Noncommutative correction term."""
    return qm * (l*(l+1)*(3*R-2)*r + R*(5*r - 8*R)) / (2 * r**5)

def V_total(r, qm, l):
    """Full potential."""
    return V_RW(r, l) + V_q(r, qm, l)

def tortoise(r, qm):
    """Tortoise coordinate r_s (Eq. 25)."""
    if r <= R + 1e-10:
        return -np.inf
    return r + R * np.log((r - R) / R) + (qm / 2) * R / (r - R)

def build_V_rs(qm, l, r_min=2.001, r_max=50.0, n_points=2000):
    """
    Build interpolation V(r_s) from the original potential.
    Returns: interpolation function, rs array, V array.
    """
    r_vals = np.linspace(r_min, r_max, n_points)
    rs_vals = np.array([tortoise(r, qm) for r in r_vals])
    V_vals = np.array([V_total(r, qm, l) for r in r_vals])
    # sort by rs
    idx = np.argsort(rs_vals)
    rs_vals = rs_vals[idx]
    V_vals = V_vals[idx]
    V_interp = interp1d(rs_vals, V_vals, kind='cubic', fill_value=0.0, bounds_error=False)
    return V_interp, rs_vals, V_vals

# ============================================================
# 2.  PROJECT 1: PLOT THE POTENTIAL (with Zeeman splitting)
# ============================================================
def plot_potential_vs_r(l=2, qm_list=[-0.2, -0.1, 0.0, 0.1, 0.2], r_max=20.0, r_min=2.01):
    """Plot V(r) vs r for multiple qm values."""
    plt.figure(figsize=(8,5))
    r_vals = np.linspace(r_min, r_max, 500)
    for qm in qm_list:
        V_vals = [V_total(r, qm, l) for r in r_vals]
        plt.plot(r_vals, V_vals, label=f'qm = {qm}')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$V$')
    plt.title(f'Regge-Wheeler potential for ℓ = {l}')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_potential_zoomed(l=2, qm_list=[-0.2, -0.1, 0.0, 0.1, 0.2], r_center=3.0, width=0.4):
    """Zoom near the peak, like Fig. 2."""
    plt.figure(figsize=(8,5))
    r_min = r_center - width
    r_max = r_center + width
    r_vals = np.linspace(r_min, r_max, 500)
    for qm in qm_list:
        V_vals = [V_total(r, qm, l) for r in r_vals]
        plt.plot(r_vals, V_vals, label=f'qm = {qm}')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$V$')
    plt.title(f'Regge-Wheeler potential (zoomed) for ℓ = {l}')
    plt.legend()
    plt.grid(True)
    plt.show()

# ============================================================
# 3.  PROJECT 2: TIME EVOLUTION (simplified 1D wave equation)
# ============================================================
def find_peak(V_interp, rs_vals):
    """Find rs where V is maximum."""
    res = minimize_scalar(lambda rs: -V_interp(rs),
                          bounds=(rs_vals[0], rs_vals[-1]),
                          method='bounded')
    return res.x, V_interp(res.x)

def time_evolution(l=2, qm=0.0, t_max=200.0, dx=0.05, cfl=0.5):
    """
    Solve the time‑dependent wave equation:
        ∂²ψ/∂t² = ∂²ψ/∂r_s² - V(r_s) ψ
    using a finite‑difference leapfrog method.
    Initial condition: Gaussian pulse.
    Returns: time array, recorded signal at a fixed observer r_s_obs.
    """
    # Build potential interpolation
    V_interp, rs_vals, _ = build_V_rs(qm, l)
    rs_min, rs_max = rs_vals[0], rs_vals[-1]
    N = int((rs_max - rs_min) / dx) + 1
    rs_grid = np.linspace(rs_min, rs_max, N)
    dx = rs_grid[1] - rs_grid[0]   # actual dx
    dt = cfl * dx
    Nt = int(t_max / dt) + 1

    # Potential on grid
    V_grid = V_interp(rs_grid)

    # Initial condition: Gaussian pulse centered at rs = -10, width 2
    psi = np.exp(-0.5 * ((rs_grid + 10) / 2)**2)
    psi_prev = psi.copy()          # ψ at t = -dt (for leapfrog)
    # For the first step, we assume ∂ψ/∂t = 0, so psi_prev = psi - dt*0 = psi

    # Observer position: peak of the potential (or fixed)
    rs_peak, _ = find_peak(V_interp, rs_vals)
    idx_obs = np.argmin(np.abs(rs_grid - rs_peak))

    signal = []
    t_vals = []

    # Leapfrog iteration
    for i in range(Nt):
        t = i * dt
        t_vals.append(t)
        signal.append(psi[idx_obs])

        # Update using finite differences
        psi_new = 2*psi - psi_prev + dt**2 * (
            (np.roll(psi, -1) - 2*psi + np.roll(psi, 1)) / dx**2 - V_grid * psi
        )
        # Simple absorbing boundaries (first‑order)
        psi_new[0] = psi_new[1] - (psi_new[2] - psi_new[1])   # crude
        psi_new[-1] = psi_new[-2] - (psi_new[-3] - psi_new[-2])

        psi_prev, psi = psi, psi_new

    # Plot signal
    plt.figure(figsize=(10,4))
    plt.plot(t_vals, np.real(signal), label='ψ at observer')
    plt.xlabel('time')
    plt.ylabel('ψ')
    plt.title(f'Time evolution (ℓ={l}, qm={qm})')
    plt.grid(True)
    plt.show()
    return t_vals, signal

# ============================================================
# 4.  PROJECT 3: QNM FREQUENCIES VIA 3rd‑ORDER WKB
# ============================================================
def derivative(f, x, dx=1e-5, n=1):
    """Numerical nth derivative (n=1..4) using central differences."""
    if n == 1:
        return (f(x+dx) - f(x-dx)) / (2*dx)
    elif n == 2:
        return (f(x+dx) - 2*f(x) + f(x-dx)) / (dx*dx)
    elif n == 3:
        return (f(x+2*dx) - 2*f(x+dx) + 2*f(x-dx) - f(x-2*dx)) / (2*dx**3)
    elif n == 4:
        return (f(x+2*dx) - 4*f(x+dx) + 6*f(x) - 4*f(x-dx) + f(x-2*dx)) / (dx**4)
    else:
        raise NotImplementedError("n>4 not implemented")

def wkb_3rd_order(V0, V2, V3, V4, n=0):
    """Return ω using 3rd‑order WKB (Iyer & Will 1987)."""
    sqrt_term = np.sqrt(-2 * V2)   # real positive
    # 2nd‑order part
    omega2 = V0 - 1j * (n + 0.5) * sqrt_term

    # Correction terms Λ2 and Λ3
    Lambda2 = (1.0 / (8.0 * V2)) * (
        (V4 / V2) * (0.25 + (n + 0.5)**2) -
        (V3**2 / V2**2) * (7.0/60.0 + (n + 0.5)**2 / 12.0)
    )
    Lambda3 = ((n + 0.5) / (16.0 * V2**2)) * (
        (V4 / V2) * (0.25 + (n + 0.5)**2) -
        (V3**2 / V2**2) * (77.0/60.0 + (n + 0.5)**2 / 12.0)
    )

    # Full 3rd‑order ω²
    omega2 += -1j * (n + 0.5) * sqrt_term * (Lambda2 + Lambda3)

    omega = np.sqrt(omega2)
    if np.real(omega) < 0:
        omega = -omega
    return omega

def compute_qnm_table(l=2, qm_list=[-0.2, -0.1, 0.0, 0.1, 0.2], n=0):
    """Print a table of QNM frequencies for given l and qm values."""
    print(f"l = {l}, n = {n}")
    print("qm      Re(ω)         -Im(ω)")
    print("----------------------------------")
    for qm in qm_list:
        V_interp, rs_vals, _ = build_V_rs(qm, l)
        rs_peak, V_peak = find_peak(V_interp, rs_vals)
        V2 = derivative(V_interp, rs_peak, dx=1e-5, n=2)
        V3 = derivative(V_interp, rs_peak, dx=1e-5, n=3)
        V4 = derivative(V_interp, rs_peak, dx=1e-5, n=4)
        omega = wkb_3rd_order(V_peak, V2, V3, V4, n)
        print(f"{qm:5.2f}   {np.real(omega):8.5f}   { -np.imag(omega):8.5f}")

# ============================================================
# 5.  MAIN: RUN ALL THREE PROJECTS (with paper‑style plots)
# ============================================================
if __name__ == "__main__":
    # Project 1: Potential plots
    print("=== Project 1: Plotting the potential ===")
    # Full range (like Fig. 1)
    plot_potential_vs_r(l=2)
    # Zoomed near peak (like Fig. 2)
    plot_potential_zoomed(l=2, r_center=3.0, width=0.4)

    # Project 3: QNM frequencies
    print("\n=== Project 3: QNM frequencies ===")
    compute_qnm_table(l=2)

  
