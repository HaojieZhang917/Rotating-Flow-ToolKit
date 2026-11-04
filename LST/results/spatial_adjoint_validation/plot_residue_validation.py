from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'residue_validation.csv').open()))
d=np.array([1e-2,3e-3,1e-3,3e-4])
fig,ax=plt.subplots(figsize=(6,4))
for prefix,label,color in [('vol','volume forcing','tab:blue'),('wall','wall forcing','tab:orange'),('full','volume+wall','tab:green')]:
    e=np.array([float(rows[-1][prefix+'_d'+str(i)]) for i in range(1,5)])
    ax.loglog(d,e,'o-',label=label,color=color)
ax.loglog(d,d/d[0]*float(rows[-1]['full_d1']),'k--',label='O(|delta|) reference')
ax.set(xlabel=r'$|\delta|$',ylabel='relative residue error',title='Spatial modal-residue validation, N=109'); ax.grid(alpha=.3,which='both'); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(root/'residue_validation.png',dpi=180); fig.savefig(root/'residue_validation.pdf'); plt.close(fig)
N=np.array([int(x['N']) for x in rows]); Q=np.array([complex(float(x['Q_real']),float(x['Q_imag'])) for x in rows])
fig,ax=plt.subplots(figsize=(6,4)); ax.plot(N,Q.real,'o-',label='Re Q'); ax.plot(N,Q.imag,'s-',label='Im Q'); ax.set(xlabel='N',ylabel='Q_discrete',title='Discrete normalization factor'); ax.grid(alpha=.3); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(root/'Q_convergence.png',dpi=180); fig.savefig(root/'Q_convergence.pdf'); plt.close(fig)
