"""K(R) and profiles for compare_turbulence_closures.jl."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results' / 'turbulent_baseflow'
rows = list(csv.DictReader((root / 'closures.csv').open()))
ST = {'ML': ('#2a78d6', '-', 'mixing length (δ_s = overshoot peak)'), 'ML99': ('#2a78d6', ':', 'mixing length (δ_s = 0.99K)'),
      'SA': ('#eb6834', '-', 'Spalart–Allmaras, confined to stator layer'),
      'SAu': ('#eb6834', '--', 'SA unconfined (ν̃ carried into rotor layer)')}
plt.rcParams.update({'axes.grid': True, 'grid.color': '#e4e3de', 'font.size': 9})
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), dpi=150)
for m, (c, ls, lab) in ST.items():
    q = sorted([r for r in rows if r['model'] == m], key=lambda r: float(r['R']))
    R = [float(r['R']) for r in q]
    ax[0].plot(R, [float(r['K']) for r in q], color=c, ls=ls, lw=1.6, label=lab)
    ax[1].plot(R, [float(r['max_nut']) for r in q], color=c, ls=ls, lw=1.6, label=lab)
ax[0].axhline(0.3147, color='#2b2b2b', lw=1.0, label='laminar Batchelor (0.3147)')
ax[0].axvline(330, color='#8a8a85', lw=0.8, ls=':'); ax[0].text(335, 0.21, 'Type-I target R=330', fontsize=7, color='#8a8a85')
ax[0].set_xlabel('R = r/δ   (Ωr²/ν = R²)'); ax[0].set_ylabel('core swirl K'); ax[0].set_title('core swirl vs radius, $Re_h$=1000'); ax[0].legend(fontsize=7, frameon=False)
ax[1].set_xlabel('R'); ax[1].set_ylabel(r'max $\nu_t/\nu$'); ax[1].set_title('eddy-viscosity level')
for m, (c, ls, lab) in ST.items():
    q = sorted([r for r in rows if r['model'] == m], key=lambda r: float(r['R']))
    ax[2].plot([float(r['R']) for r in q], [float(r['utau']) for r in q], color=c, ls=ls, lw=1.6, label=lab)
ax[2].set_xlabel('R'); ax[2].set_ylabel(r'stator $u_\tau/(\Omega r)$'); ax[2].set_title('stator friction velocity')
fig.tight_layout(); fig.savefig(root / 'closures.png'); plt.close(fig)
print('wrote', root)
