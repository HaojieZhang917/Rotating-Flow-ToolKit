from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
rows = list(csv.DictReader((ROOT / 'reconstructed_response.csv').open()))
r = np.array([float(x['r']) for x in rows])
z = np.array([float(x['z']) for x in rows])
rs = np.unique(r)
zs = np.unique(z)
nr, nz = len(rs), len(zs)

def field(name):
    values = np.array([complex(float(x[name+'_re']), float(x[name+'_im'])) for x in rows])
    return values.reshape(nr, nz)

u, v, w = field('u'), field('v'), field('w')
E = np.array([float(x['E_local']) for x in rows[::nz]])
ir0 = int(np.argmin(abs(rs)))
irp = int(np.argmax(E))

def save(fig, name):
    fig.tight_layout()
    fig.savefig(ROOT / (name + '.png'), dpi=180)
    fig.savefig(ROOT / (name + '.pdf'))
    plt.close(fig)

for ir, label in [(ir0, 'roughness center'), (irp, 'local-energy maximum')]:
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    for a, q, name in zip(ax, [u, v, w], ['u_hf', 'v_hf', 'w_hf']):
        a.plot(zs, abs(q[ir, :]), lw=2, label='magnitude')
        a.plot(zs, q[ir, :].real, '--', label='real part')
        a.plot(zs, q[ir, :].imag, ':', label='imaginary part')
        a.set(xlabel='z (wall-normal coordinate)', ylabel=name, title=f'{name}, {label}')
        a.grid(alpha=.3)
        a.legend(frameon=False, fontsize=8)
    fig.suptitle(f'Vortex-wave/roughness second-order profiles; r={rs[ir]:.5g}; R=285.36, alpha_f=0.38482, beta_f=0.07759, omega_f=-0.07759')
    save(fig, 'coupled_profiles_center' if ir == ir0 else 'coupled_profiles_energy_peak')

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(rs, E, lw=2)
ax.axvline(rs[ir0], ls='--', color='k', label='roughness center')
ax.axvline(rs[irp], ls=':', color='tab:red', label='energy maximum')
ax.set(xlabel='r (local radial coordinate)', ylabel='E_hf(r)', title='Local coupled-response energy')
ax.grid(alpha=.3)
ax.legend(frameon=False)
save(fig, 'coupled_local_energy')

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a, q, name in zip(ax, [u, v, w], ['|u_hf|', '|v_hf|', '|w_hf|']):
    im = a.pcolormesh(rs, zs, abs(q).T, shading='auto')
    a.set(xlabel='r (local radial coordinate)', ylabel='z (wall-normal coordinate)', title=name)
    fig.colorbar(im, ax=a)
fig.suptitle('Complete complex second-order response magnitude: q_hf=q_volume+q_wall')
save(fig, 'coupled_response_profiles_2d')
print('center r=', rs[ir0], 'energy peak r=', rs[irp], 'Emax=', E[irp])
