"""Figure for compute_baseflow_sensitivity.jl: two-path comparison for laminar vs turbulent-stator base flows."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results' / 'baseflow_sensitivity'
rows = list(csv.DictReader((root / 'scan.csv').open()))
TG = [('TypeI_R330_n16', 'Type-I  R=330, $n_t$=16, $\\bar\\omega$=8', 330.0),
      ('TypeII_R100_n-7', 'Type-II  R=100, $n_t$=−7, $\\bar\\omega$=8', 100.0)]
MD = [('laminar', '#2b2b2b', 'laminar Batchelor'), ('mixing_length', '#2a78d6', 'turbulent stator, mixing length'),
      ('sa', '#eb6834', 'turbulent stator, SA')]
plt.rcParams.update({'axes.grid': True, 'grid.color': '#e4e3de', 'font.size': 9})
fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=150)
for i, (t, title, R) in enumerate(TG):
    for j, nf in enumerate((-10, -9)):
        ax = axs[i, j]
        for m, c, lab in MD:
            q = sorted([r for r in rows if r['target'] == t and r['model'] == m and int(r['n_f']) == nf],
                       key=lambda r: float(r['alpha_f']))
            if not q:
                continue
            K = float(q[0]['K'])
            ax.semilogy([float(r['alpha_f']) for r in q], [float(r['hstar_delta']) for r in q], color=c, lw=1.6,
                        label=f'{lab} (K={K:.3f}, $|\\Lambda_{{dir}}|$={float(q[0]["absL_dir"]):.2e})')
        a = 30 / R
        ax.axvspan(-a, a, color='#d9d8d2', alpha=0.5, lw=0)
        ax.set_xlabel('$\\alpha_f$'); ax.set_ylabel('parity roughness height $h^*/\\delta$')
        ax.set_title(f'{title}, $n_f$={nf} (c={8/abs(nf):.3f}), forcing $e_r$', fontsize=9)
        ax.legend(fontsize=7, frameon=False)
fig.suptitle('Base-flow sensitivity of the direct vs roughness-scattering comparison ($Re_h$=1000; grey: $|\\alpha_f|R<30$)')
fig.tight_layout(); fig.savefig(root / 'baseflow_sensitivity.png'); fig.savefig(root / 'baseflow_sensitivity.pdf'); plt.close(fig)
print('wrote', root)
