"""Rotor-layer growth maps (compute_cavity_rotor_map.jl) with the gap-transmission window overlaid.

Incident stator disturbance in the long-wave limit alpha_f -> 0: k_h = |n_f|/R, n_f = -omega_bar/c, so
T(omega_bar; c) = 1/cosh(|omega_bar| Dzeta/(c R)); transmission window |omega_bar| <~ c R/Dzeta.
Dzeta = 22.8 (Re_h = 1000, distance between the G extrema of the rotor and stator layers, delta units).
"""
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results' / 'cavity_rotor_map'
DZ = 22.8
IMIN = -0.05   # physical rotor instabilities: -0.05 < Im alpha < 0 (gap in the Im-alpha histogram; spurious near-wall branch at Im ~ -0.29)
plt.rcParams.update({'axes.grid': True, 'grid.color': '#e4e3de', 'font.size': 9})
fig, axs = plt.subplots(2, 2, figsize=(13, 9), dpi=150, gridspec_kw={'width_ratios': [2.2, 1]})
summary = []
for i, R in enumerate((330, 150)):
    rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader((root / f'candidates_R{R}.csv').open())]
    rows = [r for r in rows if r['alpha_i'] > IMIN]
    best = {}
    for r in rows:
        key = (int(r['n']), int(r['wbar']))
        if key not in best or r['alpha_i'] < best[key]['alpha_i']:
            best[key] = r
    ns = sorted({k[0] for k in best}); ws = sorted({k[1] for k in best})
    Z = np.full((len(ws), len(ns)), np.nan)
    for (n, w), r in best.items():
        if r['alpha_i'] < 0:
            Z[ws.index(w), ns.index(n)] = -r['alpha_i']
    ax = axs[i, 0]
    pc = ax.pcolormesh(ns, ws, Z, shading='nearest', cmap='Blues')
    fig.colorbar(pc, ax=ax, label=r'max spatial growth $-\mathrm{Im}\,\alpha$ (rotor-layer modes)')
    for c, ls in ((0.685, ':'), (0.8, '--'), (1.0, '-')):
        wmax = c * R / DZ
        ax.axhline(wmax, color='#eb6834', ls=ls, lw=1.2, label=f'window edge $|\\bar\\omega|=cR/\\Delta\\zeta$, c={c}')
        ax.axhline(-wmax, color='#eb6834', ls=ls, lw=1.2)
    ax.set_xlabel('$n_t$ (target azimuthal number)'); ax.set_ylabel(r'$\bar\omega_t$ (rotor frame)')
    ax.set_title(f'R = {R}: unstable rotor modes and the gap-transmission window ($Re_h$=1000)')
    ax.legend(fontsize=7, frameon=False, loc='lower right')
    # frequency-wise: max growth vs omega and transmission T(omega)
    ax2 = axs[i, 1]
    gmax = [np.nanmax(Z[j]) if np.any(np.isfinite(Z[j])) else np.nan for j in range(len(ws))]
    ax2.plot(gmax, ws, color='#2a78d6', lw=1.6, label=r'max growth over $n_t$')
    ax2.set_xlabel(r'max $-\mathrm{Im}\,\alpha$'); ax2.set_ylabel(r'$\bar\omega_t$')
    ax3 = ax2.twiny()
    wq = np.linspace(min(ws), max(ws), 200)
    for c, ls in ((0.685, ':'), (0.8, '--'), (1.0, '-')):
        ax3.plot(1 / np.cosh(np.abs(wq) * DZ / (c * R)), wq, color='#eb6834', ls=ls, lw=1.2, label=f'T, c={c}')
    ax3.set_xlabel(r'gap transmission $T=1/\cosh(|\bar\omega|\Delta\zeta/(cR))$', color='#eb6834'); ax3.set_xlim(0, 1.05)
    ax3.legend(fontsize=7, frameon=False, loc='lower left'); ax2.set_title(f'R = {R}', fontsize=9)
    unst = sorted({w for (n, w), r in best.items() if r['alpha_i'] < 0})
    summary.append((R, unst, [(w, gmax[ws.index(w)]) for w in unst]))
fig.tight_layout(); fig.savefig(root / 'rotor_map_window.png'); fig.savefig(root / 'rotor_map_window.pdf'); plt.close(fig)
for R, unst, g in summary:
    print(f'R={R}: unstable omega_bar = {unst}')
    print('   max growth per omega:', ', '.join(f'{w}:{v:.4f}' for w, v in g))
    for c in (0.685, 0.8, 1.0):
        wmax = c * R / DZ
        ov = [w for w in unst if abs(w) <= wmax]
        print(f'   c={c}: window |omega|<={wmax:.1f} -> overlap {ov[0] if ov else None} .. {ov[-1] if ov else None}')
print('wrote', root)
