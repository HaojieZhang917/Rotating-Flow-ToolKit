import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out=os.path.join(os.path.dirname(__file__),"results","negative_branch_peak")
def read(name):
    with open(os.path.join(out,name),newline="") as f: return list(csv.DictReader(f))
def a(rows,key): return np.array([float(r[key]) for r in rows])
def save(fig,name):
    fig.savefig(os.path.join(out,name+".png"),dpi=220); fig.savefig(os.path.join(out,name+".pdf")); plt.close(fig)
note="alpha_f=0.38482, beta_f=0.07759, omega_f=-0.07759, signed kappa=-|kappa|, A_inc=1e-3, epsilon=1e-3"
r=read("negative_branch_refined.csv"); k=a(r,"kappa"); gf=a(r,"Gf"); gv=a(r,"Gv"); gw=a(r,"Gw"); src=a(r,"source_norm"); he=a(r,"H_eff")
imax=int(np.argmax(gf)); x=k[imax]
fig,ax=plt.subplots(figsize=(8,5)); ax.plot(k,gf,"o-",label="G_full"); ax.plot(k,gv,"--",label="G_vol"); ax.plot(k,gw,"--",label="G_wall"); ax.axvline(x,color="k",ls=":",label=f"grid peak κ={x:.3f}"); ax.set_xlabel(r"$\kappa$ (positive magnitude; solver uses $-\kappa$)"); ax.set_ylabel("normalized gain"); ax.set_title("Negative-branch refined scan"); ax.grid(alpha=.25); ax.legend(); fig.text(.02,.01,note,fontsize=7); fig.tight_layout(); save(fig,"figure1_negative_peak")

c=read("negative_branch_convergence.csv")
fig,axs=plt.subplots(1,2,figsize=(12,4))
for kk in sorted(set(float(z["kappa"]) for z in c)):
    rr=[z for z in c if float(z["kappa"])==kk]
    for N,sty in [(109,"o-"),(129,"s--"),(169,"^:")]:
        z=[q for q in rr if int(q["N"])==N]; z=sorted(z,key=lambda q:float(q["zmax"]))
        axs[0].plot([float(q["zmax"]) for q in z],[float(q["Gf"]) for q in z],sty,label=f"κ={kk:g}, N={N}")
    z=[q for q in rr if float(q["zmax"])==30.0]; z=sorted(z,key=lambda q:int(q["N"]))
    axs[1].plot([int(q["N"]) for q in z],[float(q["Gf"]) for q in z],"o-",label=f"κ={kk:g}")
axs[0].set_xlabel("zmax"); axs[1].set_xlabel("N at zmax=30"); axs[0].set_ylabel("G_full"); axs[1].set_ylabel("G_full"); axs[0].grid(alpha=.25); axs[1].grid(alpha=.25); axs[1].legend(fontsize=7)
fig.suptitle("Peak amplitude convergence",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure2_peak_convergence")

fig,axs=plt.subplots(3,1,figsize=(8,9),sharex=True)
axs[0].plot(k,gv,label="G_vol"); axs[0].plot(k,gf,label="G_full"); axs[0].set_ylabel("gain"); axs[0].legend(); axs[0].grid(alpha=.25)
axs[1].plot(k,src,label="||S_vol||"); axs[1].set_yscale("log"); axs[1].set_ylabel("source norm"); axs[1].legend(); axs[1].grid(True,which="both",alpha=.25)
axs[2].plot(k,he,label="H_eff"); axs[2].set_xlabel(r"$\kappa$"); axs[2].set_ylabel(r"$\sqrt{E_{vol}}/||S||$"); axs[2].legend(); axs[2].grid(alpha=.25)
fig.suptitle("Source strength and effective response amplification",y=.995); fig.text(.02,.005,note,fontsize=7); fig.tight_layout(); save(fig,"figure3_source_heff")

dx=np.gradient(k); dlogS=np.gradient(np.log(src),k); dlogH=np.gradient(np.log(he),k); dlogE=.5*np.gradient(np.log(gv),k)
fig,ax=plt.subplots(figsize=(8,5)); ax.plot(k,dlogS,label="d log S / dκ"); ax.plot(k,dlogH,label="d log H_eff / dκ"); ax.plot(k,dlogE,label="0.5 d log E_vol / dκ"); ax.axhline(0,color="k",lw=.7); ax.set_xlabel(r"$\kappa$"); ax.set_ylabel("logarithmic derivative"); ax.set_title("Peak mechanism derivative decomposition"); ax.grid(alpha=.25); ax.legend(); fig.text(.02,.01,note,fontsize=7); fig.tight_layout(); save(fig,"figure4_derivative_decomposition")

p=read("negative_branch_profiles.csv"); ks=sorted(set(float(z["kappa"]) for z in p)); names=[("S_u","Su_re","Su_im"),("S_v","Sv_re","Sv_im"),("S_w","Sw_re","Sw_im")]
fig,axs=plt.subplots(1,3,figsize=(14,4),sharex=True)
for ax,(nm,re,im) in zip(axs,names):
    for kk in ks:
        rr=[z for z in p if float(z["kappa"])==kk]; z=np.array([float(q["z"]) for q in rr]); y=np.hypot(a(rr,re),a(rr,im)); ax.plot(z,y,label=f"κ={kk:g}")
    ax.set_title(nm); ax.set_xlabel("z"); ax.grid(alpha=.25)
axs[0].set_ylabel("source magnitude"); axs[-1].legend(fontsize=7); fig.suptitle("Volume source spatial structure",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure5_source_profiles")

fig,axs=plt.subplots(1,3,figsize=(14,4),sharex=True)
for ax,(nm,re,im) in zip(axs,[("u","uf_re","uf_im"),("v","vf_re","vf_im"),("w","wf_re","wf_im")]):
    for kk in ks:
        rr=[z for z in p if float(z["kappa"])==kk]; z=np.array([float(q["z"]) for q in rr]); y=np.hypot(a(rr,re),a(rr,im)); ax.plot(z,y,label=f"κ={kk:g}")
    ax.set_title(nm); ax.set_xlabel("z"); ax.grid(alpha=.25)
axs[0].set_ylabel("full response magnitude"); axs[-1].legend(fontsize=7); fig.suptitle("Full second-order response profiles",y=1.02); fig.text(.02,-.02,note,fontsize=7); fig.tight_layout(); save(fig,"figure6_response_profiles")

fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
rr=sorted(r,key=lambda z:float(z["kappa"])); axs[0].plot(k,a(rr,"Rwall"),label="R_wall"); axs[0].set_yscale("log"); axs[0].set_ylabel(r"$E_{wall}/E_{vol}$"); axs[0].grid(True,which="both",alpha=.25); axs[0].legend()
axs[1].plot(k,a(rr,"rho"),label="rho"); axs[1].axhline(0,color="k",lw=.7); axs[1].set_xlabel(r"$\kappa$"); axs[1].set_ylabel(r"$\rho$"); axs[1].grid(alpha=.25); axs[1].legend()
fig.suptitle("Wall contribution and volume-wall interference",y=.995); fig.text(.02,.005,note,fontsize=7); fig.tight_layout(); save(fig,"figure7_wall_interference")
print("negative peak plots written")
