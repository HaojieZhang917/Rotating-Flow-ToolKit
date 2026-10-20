"""Read-only training-label and fixed-validation diagnostics (never uses test labels)."""
from pathlib import Path
import json
import numpy as np
import torch
from scipy.spatial import cKDTree
from scipy.interpolate import PchipInterpolator
from train_mlp import MLP

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs" / "ts_0.0" / "diagnostics"
OUT.mkdir(parents=True, exist_ok=True)
torch.set_num_threads(2)
with np.load(ROOT / "data" / "ts_0.0" / "dataset.npz", allow_pickle=False) as d:
    # Keep only the two authorized diagnostics splits in memory.
    train_idx, val_idx = d["train_idx"], d["val_idx"]
    tx, ty = d["x_raw"][train_idx], d["y_raw"][train_idx]
    vx, vy = d["x_raw"][val_idx], d["y_raw"][val_idx]
with np.load(ROOT / "data" / "ts_0.0" / "scaler.npz", allow_pickle=False) as d:
    scales = {name: d[name] for name in d.files}

def metrics(truth, prediction):
    error = prediction - truth
    return {"n": len(truth), "rmse": np.sqrt(np.mean(error**2, axis=0)).tolist(),
            "mae": np.mean(np.abs(error), axis=0).tolist(),
            "abs_error_quantiles": np.quantile(np.abs(error), [0.5, .9, .95, .99, .999, 1], axis=0).tolist()}

ckpt = torch.load(ROOT / "runs" / "ts_0.0" / "baseline" / "best_model.pt", map_location="cpu", weights_only=True)
model = MLP(ckpt["architecture"]["hidden_sizes"])
model.load_state_dict(ckpt["state_dict"])
model.eval()
with torch.no_grad():
    scaled = ((vx-scales["x_mean"])/scales["x_std"]).astype(np.float32)
    pred = np.concatenate([model(torch.from_numpy(part)).numpy() for part in np.array_split(scaled, 8)])
    pred = pred.astype(np.float64)*scales["y_std"] + scales["y_mean"]
error = pred - vy
report = {"splits_used": ["train", "validation"], "component_order": ["omega_r", "omega_i"],
          "baseline_validation": metrics(vy, pred), "groups": {}, "label_adjacent_jumps": {},
          "top_validation_errors": {}}
for j, name in enumerate(["R", "alpha", "beta"]):
    edges = np.linspace(vx[:, j].min(), vx[:, j].max(), 9)
    bins = np.clip(np.searchsorted(edges, vx[:, j], side="right")-1, 0, 7)
    groups = []
    for k in range(8):
        mask = bins == k
        groups.append({"range": [float(edges[k]), float(edges[k+1])],
                       **metrics(vy[mask], pred[mask]),
                       "squared_error_fraction": (np.sum(error[mask]**2, axis=0)/np.sum(error**2, axis=0)).tolist()})
    report["groups"][name] = groups
for j, name in enumerate(["omega_r", "omega_i"]):
    order = np.argsort(np.abs(error[:, j]))[-25:][::-1]
    report["top_validation_errors"][name] = [{"dataset_index": int(val_idx[k]), "input": vx[k].tolist(),
              "true": vy[k].tolist(), "pred": pred[k].tolist(), "error": error[k].tolist()} for k in order]
    squared = np.sort(error[:, j]**2)[::-1]
    report["baseline_validation"][name+"_sse_top_fraction"] = {
        str(f): float(squared[:max(1,int(f*len(squared)))].sum()/squared.sum()) for f in [.001,.01,.05,.1]}

# Only training labels are used for this continuity inspection.
for axis, (name, step) in enumerate(zip(["R", "alpha", "beta"], [5.0,.01,.005])):
    other = [j for j in range(3) if j != axis]
    order = np.lexsort((tx[:,axis],tx[:,other[1]],tx[:,other[0]]))
    xx, yy = tx[order], ty[order]
    gap = xx[1:,axis]-xx[:-1,axis]
    adjacent = (np.all(xx[1:,other] == xx[:-1,other], axis=1) & (gap > 0) & (gap <= step*1.001))
    starts = np.flatnonzero(adjacent)
    delta = yy[starts+1]-yy[starts]
    abs_delta = np.abs(delta)
    largest = np.argsort(np.max(abs_delta/np.array([.02,.002]),axis=1))[-25:][::-1]
    report["label_adjacent_jumps"][name] = {
        "n_adjacent_pairs": len(starts), "max_gap": step,
        "absolute_delta_quantiles": np.quantile(abs_delta,[.5,.9,.99,.999,1],axis=0).tolist(),
        "counts_delta_wr_gt_0.02": int((abs_delta[:,0]>.02).sum()),
        "counts_delta_wi_gt_0.002": int((abs_delta[:,1]>.002).sum()),
        "largest_pairs": [{"x0":xx[starts[k]].tolist(),"x1":xx[starts[k]+1].tolist(),
                           "y0":yy[starts[k]].tolist(),"y1":yy[starts[k]+1].tolist(),
                           "delta":delta[k].tolist()} for k in largest]}

# Geometric reference: train-only nearest neighbours and local weighted linear fit.
tz = (tx-scales["x_mean"])/scales["x_std"]
vz = (vx-scales["x_mean"])/scales["x_std"]
tree = cKDTree(tz)
dist, indices = tree.query(vz, k=16, workers=2)
nearest = ty[indices[:,0]]
weight = 1/np.maximum(dist[:,:8],1e-12)**2
idw = np.sum(ty[indices[:,:8]]*weight[:,:,None],axis=1)/weight.sum(axis=1)[:,None]
design = np.concatenate([np.ones((*indices.shape,1)),tz[indices]-vz[:,None,:]],axis=2)
weight = 1/np.maximum(dist,1e-12)**2
weight /= weight.max(axis=1)[:,None]
gram = np.einsum("nki,nk,nkj->nij",design,weight,design)
gram[:,np.arange(4),np.arange(4)] += 1e-10
rhs = np.einsum("nki,nk,nkc->nic",design,weight,ty[indices])
local_linear = np.linalg.solve(gram,rhs)[:,0,:]
report["train_only_interpolation"] = {"nearest": metrics(vy,nearest), "idw_8":metrics(vy,idw),
                                      "weighted_local_linear_16":metrics(vy,local_linear)}

# Single-coordinate linear interpolation with both bracketing samples in train.
coordinates = [np.unique(tx[:,j]) for j in range(3)]
train_grid = np.column_stack([np.searchsorted(coordinates[j],tx[:,j]) for j in range(3)])
valid_grid = np.column_stack([np.searchsorted(coordinates[j],vx[:,j]) for j in range(3)])
grid = np.full(tuple(len(c) for c in coordinates),-1,dtype=np.int64)
grid[tuple(train_grid.T)] = np.arange(len(tx))
axis_prediction = np.full((len(vx),3,2),np.nan)
axis_span = np.full((len(vx),3),np.nan)
for axis,name in enumerate(["R","alpha","beta"]):
    indexshape = [1,1,1]; indexshape[axis] = grid.shape[axis]
    positions = np.arange(grid.shape[axis]).reshape(indexshape)
    left = np.maximum.accumulate(np.where(grid>=0,positions,-1),axis=axis)
    right = np.flip(np.minimum.accumulate(np.flip(np.where(grid>=0,positions,grid.shape[axis]),axis=axis),axis=axis),axis=axis)
    l = left[tuple(valid_grid.T)]; r = right[tuple(valid_grid.T)]
    covered = (l>=0)&(r<grid.shape[axis])&(r>l)
    loc = np.flatnonzero(covered)
    lc = valid_grid[loc].copy(); rc = valid_grid[loc].copy()
    lc[:,axis] = l[loc]; rc[:,axis] = r[loc]
    lt,rt = grid[tuple(lc.T)],grid[tuple(rc.T)]
    fraction = (vx[loc,axis]-tx[lt,axis])/(tx[rt,axis]-tx[lt,axis])
    axis_prediction[loc,axis] = ty[lt]+fraction[:,None]*(ty[rt]-ty[lt])
    axis_span[loc,axis] = tx[rt,axis]-tx[lt,axis]
    report["train_only_interpolation"][name+"_linear"] = {
        "coverage":float(covered.mean()),**metrics(vy[loc],axis_prediction[loc,axis]),
        "span_quantiles":np.quantile(axis_span[loc,axis],[.5,.9,.99,1]).tolist()}
weight = np.where(np.isfinite(axis_span),1/(axis_span/np.array([5.,.01,.005]))**2,0)
denom = weight.sum(axis=1)
covered = denom>0
axis_mix = np.nansum(axis_prediction*weight[:,:,None],axis=1)/np.maximum(denom,1e-30)[:,None]
report["train_only_interpolation"]["coordinate_linear_mix"] = {
    "coverage":float(covered.mean()),**metrics(vy[covered],axis_mix[covered])}

# Complete the same-(R,alpha) reference at beta endpoints using two-sided
# interpolation in the interior and the nearest two train points outside.
beta_linear_full = np.full_like(vy,np.nan)
beta_pchip_full = np.full_like(vy,np.nan)
train_groups = {}
valid_groups = {}
for i,row in enumerate(tx):
    train_groups.setdefault((row[0],row[1]),[]).append(i)
for i,row in enumerate(vx):
    valid_groups.setdefault((row[0],row[1]),[]).append(i)
extrapolated = np.zeros(len(vx),dtype=bool)
for key,vlist in valid_groups.items():
    tlist = np.asarray(train_groups[key])
    tlist = tlist[np.argsort(tx[tlist,2])]
    vlist = np.asarray(vlist)
    beta = tx[tlist,2]; target = ty[tlist]
    query = vx[vlist,2]
    right = np.clip(np.searchsorted(beta,query),1,len(beta)-1)
    left = right-1
    frac = (query-beta[left])/(beta[right]-beta[left])
    linear = target[left]+frac[:,None]*(target[right]-target[left])
    cubic = PchipInterpolator(beta,target,axis=0,extrapolate=False)(query)
    outside = (query<beta[0])|(query>beta[-1])
    cubic[outside] = linear[outside]
    beta_linear_full[vlist] = linear
    beta_pchip_full[vlist] = cubic
    extrapolated[vlist] = outside
report["train_only_interpolation"]["beta_linear_full_endpoint_extrapolation"] = {
    **metrics(vy,beta_linear_full),"endpoint_extrapolation_count":int(extrapolated.sum())}
report["train_only_interpolation"]["beta_pchip_full_endpoint_linear_extrapolation"] = {
    **metrics(vy,beta_pchip_full),"endpoint_extrapolation_count":int(extrapolated.sum())}
for j,name in enumerate(["omega_r","omega_i"]):
    good = covered
    component_error = axis_mix[good,j]-vy[good,j]
    order = np.argsort(np.abs(component_error))[-15:][::-1]
    loc = np.flatnonzero(good)[order]
    report["train_only_interpolation"]["coordinate_linear_mix"][name+"_largest_errors"] = [
        {"input":vx[k].tolist(),"true":vy[k].tolist(),"pred":axis_mix[k].tolist(),
         "axis_predictions":np.nan_to_num(axis_prediction[k],nan=0).tolist()} for k in loc]

np.savez_compressed(OUT/"validation_diagnostics.npz",dataset_indices=val_idx,x=vx,y=vy,
                    baseline_prediction=pred,nearest_prediction=nearest,idw_prediction=idw,
                    local_linear_prediction=local_linear,coordinate_prediction=axis_prediction,
                    coordinate_mix_prediction=axis_mix,coordinate_covered=covered,
                    beta_linear_full=beta_linear_full,beta_pchip_full=beta_pchip_full,
                    beta_endpoint_extrapolation=extrapolated)
(OUT/"validation_diagnostics.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
print(json.dumps({"baseline":report["baseline_validation"],
                  "interpolation":{name:{k:v for k,v in value.items() if k in ["n","rmse","mae","coverage"]}
                                   for name,value in report["train_only_interpolation"].items()},
                  "jumps":{name:{k:v for k,v in value.items() if k != "largest_pairs"}
                           for name,value in report["label_adjacent_jumps"].items()}},indent=2))
