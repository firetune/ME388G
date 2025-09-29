"""
Uranium slab spatial profiles (no scattering, mono-directional beam)
- For each neutron energy: plot
  (1) Fission neutron source and specific heat (same axes)
  (2) 239Pu production rate (from 238U capture)
  (3) Total absorption (loss) rate
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# -------------------- constants & inputs --------------------
phi0 = 1e14               # incident mono-directional flux [n/(cm^2·s)]
rho_g_cm3 = 19.1          # uranium density [g/cm^3]
rho_kg_cm3 = rho_g_cm3/1e3  # [kg/cm^3]
barn = 1e-24              # [cm^2]
NA = 6.022e23             # [1/mol]
nu_bar = 2.4              # neutrons per fission
E_fission_J = 200e6 * 1.602e-19  # 200 MeV -> J

# mixture: 5 wt% 235U, 95 wt% 238U
wf_235, wf_238 = 0.05, 0.95
M_235, M_238 = 235.043928117, 238.050786936  # g/mol

# number densities [atoms/cm^3]
N235 = (rho_g_cm3 * wf_235 / M_235) * NA
N238 = (rho_g_cm3 * wf_238 / M_238) * NA

# microscopic cross sections (barns) from your table
xs = {
    "1 MeV":   {"sa235": 1.2, "sf235": 1.0, "sa238": 0.3,   "sf238": 0.27},
    "6.7 eV":  {"sa235": 10., "sf235": 6.,  "sa238": 1500., "sf238": 0.2},
    "6.5 eV":  {"sa235": 10., "sf235": 6.,  "sa238": 10.,   "sf238": 0.0},
    "0.025 eV":{"sa235": 680.,"sf235": 585.,"sa238": 2.7,   "sf238": 0.0},
}

# spatial grid (1 cm slab, beam normal to surface)
x = np.linspace(0.0, 1.0, 1000)  # [cm]

# where to save plots
OUT = "plots"
os.makedirs(OUT, exist_ok=True)

def maybe_set_log_y(ax, *ys, ratio_threshold=100.0):
    """Turn on log-y if dynamic range across the given arrays is huge."""
    vals = []
    for y in ys:
        y = np.asarray(y)
        vals.append(y[(y > 0) & np.isfinite(y)])
    vals = np.concatenate(vals) if len(vals) else np.array([])
    if vals.size:
        r = np.max(vals) / max(np.min(vals), 1e-300)
        if r > ratio_threshold:
            ax.set_yscale("log")

# -------------------- loop energies & plot --------------------
for label, d in xs.items():
    # macroscopic cross sections [cm^-1]
    Sigma_f = (N235*d["sf235"] + N238*d["sf238"]) * barn
    Sigma_t = (N235*(d["sa235"] + d["sf235"]) + N238*(d["sa238"] + d["sf238"])) * barn
    Sigma_a238 = N238*d["sa238"] * barn
    Sigma_a = (N235*d["sa235"] + N238*d["sa238"]) * barn  # total absorption

    # uncollided, mono-directional flux profile (no scattering)
    phi_x = phi0 * np.exp(-Sigma_t * x)

    # reaction rates / sources
    R_fiss = Sigma_f * phi_x               # [fissions/cm^3·s]
    q_neut = nu_bar * R_fiss               # fission neutron source [n/cm^3·s]
    q_spec = (E_fission_J * R_fiss) / rho_kg_cm3  # specific power [W/kg]
    R_pu   = Sigma_a238 * phi_x            # 239Pu production [atoms/cm^3·s]
    R_loss = Sigma_a * phi_x               # total absorption [reactions/cm^3·s]
    
    MFP = 1.0 / Sigma_t if Sigma_t > 0 else np.inf  # mean free path [cm]
    
    #Print values of Sigma_t and Sigma_f and Sigma_a238 for each energy and sum of Sigma_a and Sigma_f
    print(f"Energy: {label}")
    print(f"  Sigma_t = {Sigma_t:.4e} cm^-1")
    print(f"  Sigma_f = {Sigma_f:.4e} cm^-1")
    print(f"  Sigma_a238 = {Sigma_a238:.4e} cm^-1")
    print(f"  Sigma_a = {Sigma_a:.4e} cm^-1")
    print(f"  Sigma_a + Sigma_f = {Sigma_a + Sigma_f:.4e} cm^-1")
    print()

    # (1) fission neutron source + specific heat generation
    fig = plt.figure()
    plt.plot(x, q_neut, label="Fission neutron source (n/cm$^3\\cdot$s)")
    plt.plot(x, q_spec, label="Specific heat from fission (W/kg)")
    plt.xlabel("Depth x (cm)")
    plt.ylabel("Value")
    plt.title(f"{label}: Fission Source & Specific Power vs Depth MFP={MFP:.2g} cm")
    maybe_set_log_y(plt.gca(), q_neut, q_spec)
    plt.legend()
    plt.grid(True, which="both", alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"{label.replace(' ','_').replace('.','p')}_fission_and_power.png"), dpi=160)
    plt.show()

    # (2) Pu-239 production
    fig = plt.figure()
    plt.plot(x, R_pu, label=r"$^{239}$Pu production rate (atoms/cm$^3\cdot$s)")
    plt.xlabel("Depth x (cm)")
    plt.ylabel("Rate")
    plt.title(f"{label}: Local $^{{239}}$Pu Production vs Depth")
    maybe_set_log_y(plt.gca(), R_pu)
    plt.legend()
    plt.grid(True, which="both", alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"{label.replace(' ','_').replace('.','p')}_pu_rate.png"), dpi=160)
    plt.show()

    # (3) total loss (absorption)
    fig = plt.figure()
    plt.plot(x, R_loss, label="Total absorption (reactions/cm$^3\\cdot$s)")
    plt.xlabel("Depth x (cm)")
    plt.ylabel("Rate")
    plt.title(f"{label}: Local Absorption (Loss) vs Depth MFP={MFP:.2g} cm")
    maybe_set_log_y(plt.gca(), R_loss)
    plt.legend()
    plt.grid(True, which="both", alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"{label.replace(' ','_').replace('.','p')}_loss_rate.png"), dpi=160)
    plt.show()
