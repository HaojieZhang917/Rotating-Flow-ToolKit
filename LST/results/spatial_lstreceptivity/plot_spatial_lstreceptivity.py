from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'spatial_lstr_eigenvalues.csv').open()))
fig,ax=plt.subplots(1,2,figsize=(10,4))
for case,color in [('target','tab:red'),('reference','tab:blue')]:
    a=[x for x in rows if x['case']==case and x['alpha_real']!='NaN']
    j=0 if case=='target' else 1
    ax[j].scatter([float(x['alpha_real']) for x in a],[float(x['alpha_imag']) for x in a],label=case,color=color)
    ax[j].axhline(0,color='k',lw=.6); ax[j].set(xlabel=r'Re($\alpha$)',ylabel=r'Im($\alpha$)',title=case); ax[j].grid(alpha=.3); ax[j].legend(frameon=False)
fig.suptitle('Spatial LST spectrum, R=285.36, beta=0.07759'); fig.tight_layout(); fig.savefig(ROOT/'spatial_alpha_spectrum.png',dpi=180); fig.savefig(ROOT/'spatial_alpha_spectrum.pdf'); plt.close(fig)

p=list(csv.DictReader((ROOT/'spatial_mode_profiles.csv').open()))
z=np.array([float(x['z']) for x in p])
fig,ax=plt.subplots(1,2,figsize=(10,4))
for q,label in [('u','u'),('v','v'),('w','w'),('p','p')]:
    f=np.array([complex(float(x[q+'_re']),float(x[q+'_im'])) for x in p]); ax[0].plot(z,abs(f),label=label)
for q,label in [('au','adjoint u'),('av','adjoint v'),('aw','adjoint w'),('ap','adjoint p')]:
    f=np.array([complex(float(x[q+'_re']),float(x[q+'_im'])) for x in p]); ax[1].plot(z,abs(f),label=label)
for a,title in zip(ax,['direct candidate mode','discrete transpose adjoint mode']):
    a.set(xlabel='z',ylabel='magnitude',title=title,xlim=(0,20)); a.grid(alpha=.3); a.legend(frameon=False,fontsize=8)
fig.suptitle('von Karman spatial LST profiles; reference omega=0, beta=0.07759, R=285.36'); fig.tight_layout(); fig.savefig(ROOT/'spatial_mode_profiles.png',dpi=180); fig.savefig(ROOT/'spatial_mode_profiles.pdf'); plt.close(fig)
print('plots written')
