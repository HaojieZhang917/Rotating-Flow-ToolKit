"""Profiles of the turbulent-stator base flow (validate_turbulent_baseflow.jl)."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results' / 'turbulent_baseflow'
rows = list(csv.DictReader((root / 'profiles.csv').open()))
cases = [('laminar', '#2b2b2b', '-', 'laminar'), ('R200', '#2a78d6', '-', 'R=200'), ('R330', '#eb6834', '-', 'R=330'),
         ('R450', '#1baf7a', '--', 'R=450'), ('R600', '#eda100', ':', 'R=600')]
plt.rcParams.update({'axes.grid': True, 'grid.color': '#e4e3de', 'font.size': 9})
fig, axs = plt.subplots(1, 5, figsize=(16, 4.6), dpi=150)
for key, c, ls, lab in cases:
    q = [r for r in rows if r['case'] == key]
    z = [float(r['z']) for r in q]
    for ax, var in zip(axs[:4], ('G', 'F', 'H', 'nut')):
        ax.plot([float(r[var]) for r in q], z, color=c, ls=ls, lw=1.5, label=lab)
    axs[4].plot([float(r['G']) for r in q], z, color=c, ls=ls, lw=1.5, label=lab)
for ax, t in zip(axs, ('G (lab frame)', 'F', 'H', r'$\nu_t/\nu$', 'G, stator side (zoom)')):
    ax.set_title(t); ax.set_ylabel('z (rotor 0, stator 1)')
axs[4].set_ylim(0.7, 1.0); axs[4].set_xlim(-0.02, 0.45)
axs[0].legend(fontsize=8, frameon=False)
fig.suptitle('$Re_h$=1000: base flow with turbulent stator layer (mixing length, local similarity) vs laminar Batchelor')
fig.tight_layout(); fig.savefig(root / 'turbulent_baseflow.png'); plt.close(fig)
print('wrote', root)
