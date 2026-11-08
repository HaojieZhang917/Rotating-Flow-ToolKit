import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent/'results'/'stator_forcing_demo'
meta = dict(l.strip().split('=') for l in (root/'meta.txt').open())
zre, zse, zs = float(meta['z_rotor_edge']), float(meta['z_stator_edge']), float(meta['zs'])
f = lambda r, k: float(r[k])
bf = list(csv.DictReader((root/'baseflow.csv').open()))
pr = list(csv.DictReader((root/'profiles.csv').open()))
pe = list(csv.DictReader((root/'penetration.csv').open()))

fig = plt.figure(figsize=(15, 8.2), dpi=150)
gs = fig.add_gridspec(2, 4)
# (a) base flow + forcing location
a = fig.add_subplot(gs[:, 0])
z = [f(r, 'z') for r in bf]
a.plot([f(r, 'F')*3 for r in bf], z, label='3F'); a.plot([f(r, 'G') for r in bf], z, label='G')
a.plot([f(r, 'H')/20 for r in bf], z, label='H/20')
zz = np.linspace(0, 1, 400); a.fill_betweenx(zz, 0, 0.9*np.exp(-(zz - zs)**2/(2*float(meta['sigma'])**2)), color='tab:red', alpha=0.25, label='forcing envelope')
for y, lab in ((zre, 'rotor-layer edge'), (zse, 'stator-layer edge')):
    a.axhline(y, color='k', ls=':', lw=0.8); a.text(0.55, y + 0.01, lab, fontsize=7)
a.set_xlabel('base flow'); a.set_ylabel('z (rotor z=0, stator z=1)'); a.set_title(r'$Re_h=1000$, base flow'); a.legend(fontsize=7, frameon=False)
# (b-d) response profiles for the three forcing directions
for j, d in enumerate(('r', 'theta', 'z')):
    a = fig.add_subplot(gs[0, j+1])
    for af, c in ((0.02, 'tab:blue'), (0.1, 'tab:orange'), (0.3, 'tab:green')):
        q = [r for r in pr if r['dir'] == d and f(r, 'alpha_f') == af]
        zp = [f(r, 'z') for r in q]
        a.semilogx([f(r, 'u') for r in q], zp, '-', color=c, lw=1.2, label=rf'$\alpha_f={af}$  |u|')
        a.semilogx([f(r, 'v') for r in q], zp, '--', color=c, lw=1.0, label='|v|')
        a.semilogx([f(r, 'w') for r in q], zp, ':', color=c, lw=1.0, label='|w|')
    for y in (zre, zse): a.axhline(y, color='k', ls=':', lw=0.8)
    a.axhline(0.5, color='gray', lw=0.6)
    a.set_xlim(1e-4, 1e3); a.set_title(rf'forcing $e_{{{d}}}$, $R=330$, $n_f=-10$, $c=0.85$', fontsize=9)
    a.set_xlabel('velocity amplitude (normalised at z=0.5)'); a.set_ylabel('z')
    if j == 0: a.legend(fontsize=6, frameon=False, ncol=1)
# (e-g) penetration: rotor-edge / core amplitude versus alpha_f
for j, d in enumerate(('r', 'theta', 'z')):
    a = fig.add_subplot(gs[1, j+1])
    for nf, col in ((-20, 'tab:purple'), (-10, 'tab:blue'), (0, 'tab:green'), (10, 'tab:orange')):
        for cc, ls in ((0.7, ':'), (0.85, '-'), (1.0, '--')):
            q = sorted([r for r in pe if r['dir'] == d and f(r, 'n_f') == nf and f(r, 'c') == cc], key=lambda r: f(r, 'alpha_f'))
            a.semilogy([f(r, 'alpha_f') for r in q], [f(r, 'edge_over_core') for r in q], ls, color=col, lw=1.1,
                       label=(rf'$n_f={nf}$' if cc == 0.85 else None))
    a.axvspan(0.2, 0.5, color='tab:red', alpha=0.07, label='rotor Type-I/II $\\alpha_r$ range')
    a.set_xlabel(r'$\alpha_f$'); a.set_ylabel(r'$|u|_{\rm rotor\ edge}/|u|_{z=0.5}$')
    a.set_title(rf'penetration, forcing $e_{{{d}}}$ (line style: $c$=0.7 dotted, 0.85 solid, 1 dashed)', fontsize=8)
    a.grid(True, which='both', alpha=.3)
    if j == 0: a.legend(fontsize=7, frameon=False)
fig.suptitle('Stator-side disturbance model: localised Reynolds-stress forcing in the stator layer, response across the gap')
fig.tight_layout()
fig.savefig(root/'stator_forcing_demo.png'); fig.savefig(root/'stator_forcing_demo.pdf'); plt.close(fig)
print('wrote', root)
