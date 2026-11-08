"""Figure for compute_cavity_scattering_scan.jl: direct vs roughness-scattering excitation in the cavity."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent / 'results' / 'cavity_scattering_scan'
rows = list(csv.DictReader((root / 'scan.csv').open()))
for r in rows:
    for k, v in r.items():
        try:
            r[k] = float(v)
        except ValueError:
            pass
TARGETS = [('TypeI_R330_n16', 'Type-I  R=330, $n_t$=16, $\\bar\\omega$=8'),
           ('TypeII_R100_n-7', 'Type-II  R=100, $n_t$=−7, $\\bar\\omega$=8')]
NF = [-8, -9, -10, -11]
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']          # fixed categorical order (n_f)
DIRS = [('r', '#2a78d6', '-'), ('theta', '#eb6834', '--'), ('z', '#1baf7a', ':')]
INK, MUTED = '#2b2b2b', '#8a8a85'
plt.rcParams.update({'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.grid': True, 'grid.color': '#e4e3de', 'grid.linewidth': 0.6, 'font.size': 9})


def sel(t, nf, d='r'):
    s = sorted([r for r in rows if r['target'] == t and r['n_f'] == nf and r['dir'] == d], key=lambda r: r['alpha_f'])
    return s, np.array([r['alpha_f'] for r in s])


def shade(ax, R):
    a = 30 / R
    ax.axvspan(-a, a, color='#d9d8d2', alpha=0.5, lw=0)


def dirty(s):
    return np.array([(r['d_f'] < 0.02) or (r['d_h'] < 0.02) for r in s])


fig, axs = plt.subplots(2, 4, figsize=(16, 7.6), dpi=150)
for i, (t, title) in enumerate(TARGETS):
    R = [r['R'] for r in rows if r['target'] == t][0]
    # (a) parity height
    ax = axs[i, 0]
    for nf, c in zip(NF, COL):
        s, af = sel(t, nf); h = np.array([r['hstar_delta'] for r in s]); bad = dirty(s)
        ax.semilogy(af, h, color=c, lw=1.6, label=f"$n_f$={nf}  (c={s[0]['c']:.3f})")
        ax.plot(af[bad], h[bad], 'x', color=c, ms=4, mew=1.0)
    shade(ax, R)
    ax.set_xlabel('$\\alpha_f$'); ax.set_ylabel('parity roughness height $h^*/\\delta$')
    ax.set_title(f'{title}\n$h^*=|\\Lambda_{{dir}}|/|\\Lambda_{{hf}}|$, forcing $e_r$ (× = near resonance)', fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    # (b) incident amplitude at the rotor-layer edge
    ax = axs[i, 1]
    for nf, c in zip(NF, COL):
        s, af = sel(t, nf)
        ax.semilogy(af, [r['qf_edge'] for r in s], color=c, lw=1.6, label=f'$n_f$={nf}')
    shade(ax, R)
    ax.set_xlabel('$\\alpha_f$'); ax.set_ylabel('$|q_f|$ at rotor-layer edge (unit forcing)')
    ax.set_title('incident field reaching the rotor', fontsize=9); ax.legend(fontsize=7, frameon=False)
    # (c) volume vs wall contribution, n_f = -10
    ax = axs[i, 2]
    s, af = sel(t, -10)
    for key, c, ls, lab in (('absL_hf', INK, '-', 'total $|\\Lambda_{hf}|$'), ('absL_vol', '#2a78d6', '--', 'volume'),
                            ('absL_wall', '#eb6834', ':', 'wall')):
        ax.semilogy(af, [r[key] for r in s], color=c, ls=ls, lw=1.6, label=lab)
    shade(ax, R)
    ax.set_xlabel('$\\alpha_f$'); ax.set_ylabel('$|\\Lambda|$ (unit forcing, unit $h$ in gap units)')
    ax.set_title('scattering coefficient split, $n_f$=−10, $e_r$', fontsize=9); ax.legend(fontsize=7, frameon=False)
    # (d) forcing direction, n_f = -10
    ax = axs[i, 3]
    for d, c, ls in DIRS:
        s, af = sel(t, -10, d)
        ax.semilogy(af, [r['hstar_delta'] for r in s], color=c, ls=ls, lw=1.6, label=f'$e_{{{d}}}$')
    shade(ax, R)
    ax.set_xlabel('$\\alpha_f$'); ax.set_ylabel('$h^*/\\delta$')
    ax.set_title('forcing direction, $n_f$=−10', fontsize=9); ax.legend(fontsize=7, frameon=False)
fig.suptitle('Rotor–stator cavity, $Re_h$=1000: roughness-scattering vs direct excitation by stator-layer forcing (flat forcing spectrum); grey band: $|\\alpha_f|R<30$, local-parallel approximation questionable',
             color=INK)
fig.tight_layout()
fig.savefig(root / 'cavity_scattering_scan.png'); fig.savefig(root / 'cavity_scattering_scan.pdf'); plt.close(fig)
print('wrote', root)
