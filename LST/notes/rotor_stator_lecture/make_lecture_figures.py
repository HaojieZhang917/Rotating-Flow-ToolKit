"""Didactic figures for the lecture notes (f01_setup.png, f10_filter_theory.png).
Data are read from LST/results (laminar Batchelor base flow at Re_h = 1000, near-inertial transmission test)."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

S = Path(__file__).resolve().parents[2] / 'results'
OUT = Path(__file__).resolve().parent / 'figs'
INK, MUTED, BLUE, ORANGE, GREEN = '#2b2b2b', '#8a8a85', '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'axes.grid': True, 'grid.color': '#e4e3de', 'font.size': 9})

prof = [r for r in csv.DictReader((S / 'turbulent_baseflow' / 'profiles.csv').open()) if r['case'] == 'laminar']
z = np.array([float(r['z']) for r in prof]); F = np.array([float(r['F']) for r in prof])
G = np.array([float(r['G']) for r in prof]); H = np.array([float(r['H']) for r in prof])
o = np.argsort(z); z, F, G, H = z[o], F[o], G[o], H[o]
zr, zs, K = 0.124, 0.845, 0.3147

# ---------------------------------------------------------------- f01: set-up
fig, ax = plt.subplots(1, 3, figsize=(15, 5), dpi=150, gridspec_kw={'width_ratios': [1.25, 1, 1]})
a = ax[0]; a.set_xlim(0, 10); a.set_ylim(-1.2, 11.5); a.axis('off')
a.add_patch(Rectangle((0.3, -0.6), 9.4, 0.6, color='#9db8d9'))                 # rotor
a.add_patch(Rectangle((0.3, 10.0), 9.4, 0.6, color='#c9c7bf'))                 # stator
a.text(5, -1.05, 'rotor (rotates with Ω), laminar Ekman-type layer, roughness h(r,θ)', ha='center', fontsize=8)
a.text(5, 11.0, 'stator (at rest), turbulent Bödewadt layer', ha='center', fontsize=8)
a.add_patch(Rectangle((0.3, 0.0), 9.4, 1.2, color=BLUE, alpha=0.12)); a.text(9.6, 0.55, 'rotor layer ~3δ', ha='right', fontsize=7)
a.add_patch(Rectangle((0.3, 8.4), 9.4, 1.6, color=ORANGE, alpha=0.12)); a.text(9.6, 9.1, 'stator layer', ha='right', fontsize=7)
a.text(5, 4.8, 'inviscid core: solid-body rotation KΩ\n(K ≈ 0.315), weak axial flow toward the rotor', ha='center', fontsize=8, color=INK)
xs = np.linspace(0.6, 9.4, 300); a.plot(xs, 9.2 + 0.25*np.sin(3.1*xs) + 0.12*np.sin(11*xs), color=ORANGE, lw=1)
for x0 in (2.0, 5.0, 8.0):                                                      # long waves penetrate
    a.add_patch(FancyArrowPatch((x0, 8.2), (x0, 1.5), arrowstyle='-|>', mutation_scale=12, color=GREEN, lw=1.6))
a.text(2.25, 6.6, 'long waves\n(λ ≳ 4.5h)\npass: T≈1/cosh(k_hΔζ)', fontsize=7, color=GREEN)
a.add_patch(FancyArrowPatch((6.4, 8.2), (6.4, 5.8), arrowstyle='-|>', mutation_scale=12, color=MUTED, lw=1.2, ls='--'))
a.text(6.55, 6.0, 'short waves\ndecay e^{-k_hΔζ}', fontsize=7, color=MUTED)
a.annotate('', xy=(9.95, 10.0), xytext=(9.95, 0.0), arrowprops=dict(arrowstyle='<->', color=INK, lw=0.8))
a.text(9.98, 5.0, 'h = √Re_h δ', rotation=90, va='center', fontsize=7)
a.set_title('(a) rotor–stator cavity and the physical chain', fontsize=9)
b = ax[1]
b.plot(G, z, color=ORANGE, lw=1.6, label='G (lab frame)'); b.plot(3*F, z, color=BLUE, lw=1.6, label='3F'); b.plot(H/20, z, color=GREEN, lw=1.6, label='H/20')
for y, lab in ((zr, 'z_r = 0.124 (first G minimum)'), (zs, 'z_s = 0.845 (G overshoot)')):
    b.axhline(y, color=MUTED, lw=0.8, ls=':'); b.text(-0.48, y + 0.012, lab, fontsize=7, color=MUTED)
b.axhspan(zr, zs, color='#efeee9', zorder=0); b.text(0.55, 0.47, 'core', fontsize=8, color=MUTED)
b.set_xlabel('base flow (Re_h = 1000)'); b.set_ylabel('z = z*/h  (rotor 0, stator 1)'); b.legend(fontsize=7, frameon=False, loc='center right')
b.set_title('(b) laminar Batchelor base flow', fontsize=9)
c_ = ax[2]
c_.plot(1 - G, z, color=INK, lw=1.6, label='1 − G(z): local backward speed in rotor frame')
for cc, col in ((0.70, ORANGE), (0.75, GREEN), (0.85, BLUE), (1.0, MUTED)):
    c_.axvline(cc, color=col, lw=1, ls='--', label=f'c = {cc}')
    zc = [0.5*(z[i] + z[i+1]) for i in range(len(z)-1) if np.sign(1 - G[i] - cc) != np.sign(1 - G[i+1] - cc)]
    c_.plot([cc]*len(zc), zc, 'o', color=col, ms=5)
c_.axhspan(zr, zs, color='#efeee9', zorder=0)
c_.set_xlabel('speed (units of Ωr)'); c_.set_ylabel('z'); c_.set_xlim(0, 1.05)
c_.set_title('(c) critical levels σ(z)=0  ⇔  c = 1 − G(z)  (dots)', fontsize=9); c_.legend(fontsize=7, frameon=False, loc='lower left')
fig.tight_layout(); fig.savefig(OUT / 'f01_setup.png'); plt.close(fig)

# ---------------------------------------------------------------- f10: filter theory
ni = [{k: float(v) for k, v in r.items()} for r in csv.DictReader((S / 'gap_filter' / 'near_inertial.csv').open())]
beta = -10/330
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), dpi=150)
x = np.linspace(1.003, 6, 400)
ax[0].plot(x, x/np.sqrt(x**2 - 1), color=INK, lw=1.6, label=r'inviscid  $g=m/k_h=x/\sqrt{x^2-1}$')
for al, mk in ((0.03, 'o'), (0.06, 's'), (0.12, '^')):
    q = [r for r in ni if abs(r['alpha'] - al) < 1e-9]
    kh = np.hypot(al, beta)
    ax[0].plot([r['x'] for r in q], [r['me_r']/kh for r in q], mk, color=ORANGE, ms=6, mfc='none',
               label=f'exact core root (viscosity, axial flow, curvature), α={al}')
ax[0].axhline(1, color=MUTED, lw=0.8)
ax[0].set_xlabel(r'$\sigma/f$  (core-frame frequency / inertial frequency)'); ax[0].set_ylabel(r'$m/k_h$')
ax[0].set_ylim(0, 6); ax[0].legend(fontsize=7, frameon=False)
ax[0].set_title('(a) vertical decay rate relative to potential decay', fontsize=9)
kd = np.linspace(0, 12, 400)
ax[1].semilogy(kd, 1/np.cosh(kd), color=INK, lw=1.6, label=r'$T=1/\cosh(k_h\Delta\zeta)$')
ax[1].semilogy(kd, 2*np.exp(-kd), color=MUTED, lw=1, ls='--', label=r'$2e^{-k_h\Delta\zeta}$')
pts = [(0.83, 'incident long wave\nR=330, n_f=−10, α_f=0.02', GREEN), (1.9, 'incident for Type-II\nR=150, ω̄=10, c=0.8', GREEN),
       (6.57, 'direct path Type-II\nR=100 (k_h=0.288)', ORANGE), (10.3, 'direct path Type-I\nR=330 (k_h=0.451)', ORANGE)]
for xv, lab, col in pts:
    ax[1].plot(xv, 1/np.cosh(xv), 'o', color=col, ms=7)
    ax[1].annotate(lab, (xv, 1/np.cosh(xv)), textcoords='offset points', xytext=(6, 6), fontsize=7, color=col)
ax[1].axvspan(0, 1, color=GREEN, alpha=0.08); ax[1].text(0.1, 2e-5, 'transmission\nwindow', fontsize=7, color=GREEN)
ax[1].set_xlabel(r'$k_h\Delta\zeta$  ($\Delta\zeta\approx0.72\sqrt{Re_h}$; 22.8 at $Re_h$=1000)'); ax[1].set_ylabel('T')
ax[1].set_ylim(1e-5, 2); ax[1].legend(fontsize=7, frameon=False, loc='upper right')
ax[1].set_title('(b) the gap filter: one law from long waves to the direct path', fontsize=9)
fig.tight_layout(); fig.savefig(OUT / 'f10_filter_theory.png'); plt.close(fig)
print('wrote', OUT)
