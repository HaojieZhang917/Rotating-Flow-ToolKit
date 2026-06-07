#!/usr/bin/env python3
"""
================================================================================
  动坐标系 von Kármán 旋转圆盘基本流求解器  (Boussinesq 近似, solve_bvp)
================================================================================

  方程 (7 变量 ODE 系统):
    y[0]=H, y[1]=F', y[2]=F, y[3]=G', y[4]=G, y[5]=T', y[6]=T

    H'  = -2F
    F'' = F² + H F' - (G-1)² + (T-1)
    G'' = 2FG + H G' - 2F
    T'' = Pr * H * T'

  边界条件:
    z=0:  H=0, F=0, G=0, T=Tw
    z→∞: F=0, G=1, T=1

  适用范围:
    Boussinesq 近似, Tw ∈ [0.25, 1.049]
    冷壁 (Tw<1): 大步长延拓即可
    热壁 (1<Tw≤1.049): 需小步长延拓
    Tw ≥ 1.05:  相似性解不存在 (需全耦合可压缩模型)

  用法:
    python Bone.py              # 使用文件顶部的 Tw 值
    python Bone.py 1.03         # 命令行指定 Tw (需从 Tw=1 延拓)

  输出:
    baseflow_Tw{Tw}.dat  — 基本流剖面数据
    Bone_profile_Tw{Tw}.png — 剖面图
================================================================================
"""

import sys
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ═══════════════════════════════════════════════════════════════
#  用户参数 — 修改这里
# ═══════════════════════════════════════════════════════════════
TARGET_TW  = 1.04       # 目标壁温 (0.25 ~ 1.049)
Pr         = 0.72        # 普朗特数 (空气 ≈ 0.72)
ZMAX       = 20.0        # 计算域高度
N_INITIAL  = 500         # 初始网格点数
TOL        = 1e-8        # 求解容差
MAX_NODES  = 200000      # 最大网格点数
SAVE_PLOT  = True        # 是否保存图片

# ═══════════════════════════════════════════════════════════════
#  从命令行覆盖 Tw (可选)
# ═══════════════════════════════════════════════════════════════
if len(sys.argv) > 1:
    TARGET_TW = float(sys.argv[1])
    print(f"[命令行] 覆盖 Tw = {TARGET_TW}")

# ═══════════════════════════════════════════════════════════════
#  ODE 系统 (动坐标系, Boussinesq 近似)
# ═══════════════════════════════════════════════════════════════
def vonKarman_ODE(z, y):
    """
    动坐标系可压缩 von Kármán 方程
    y[0]=H, y[1]=F', y[2]=F, y[3]=G', y[4]=G, y[5]=T', y[6]=T
    """
    return np.array([
        -2.0 * y[2],                                      # H'  = -2F
        y[2]**2 + y[0]*y[1] - (y[4]-1.0)**2 + y[6] - 1.0, # F'' = F²+HF'-(G-1)²+(T-1)
        y[1],                                              # F'  = F'
        2*y[2]*y[4] + y[0]*y[3] - 2*y[2],                 # G'' = 2FG+HG'-2F
        y[3],                                              # G'  = G'
        Pr * y[0] * y[5],                                  # T'' = Pr·H·T'
        y[5],                                              # T'  = T'
    ])

def vonKarman_bc(Tw):
    """返回给定 Tw 的边界条件函数"""
    def bc(ya, yb):
        return np.array([
            ya[0],          # H(0)  = 0
            ya[2],          # F(0)  = 0
            ya[4],          # G(0)  = 0
            ya[6] - Tw,     # T(0)  = Tw
            yb[2],          # F(∞)  = 0
            yb[4] - 1.0,    # G(∞)  = 1
            yb[6] - 1.0,    # T(∞)  = 1
        ])
    return bc

# ═══════════════════════════════════════════════════════════════
#  初始猜测 (基于经典 von Kármán 不可压解)
# ═══════════════════════════════════════════════════════════════
def initial_guess(z):
    """
    Tw=1 等温解的初始猜测。
    对于其他 Tw，将通过延拓自动调整。
    """
    g = np.zeros((7, len(z)))
    g[0] = -0.8845 * (1.0 - np.exp(-0.8 * z))           # H:  0 → -0.8845
    g[1] =  0.5102 * (1.0 - 0.8 * z) * np.exp(-0.8 * z) # F': 0.51 → 0
    g[2] =  0.5102 * z * np.exp(-0.8 * z)                # F:  0 → 0
    g[3] =  0.6159 * np.exp(-0.8 * z)                    # G': 0.616 → 0
    g[4] =  1.0 - np.exp(-0.8 * z)                       # G:  0 → 1
    g[5] =  0.0                                           # T': 0
    g[6] =  1.0                                           # T:  1
    return g

# ═══════════════════════════════════════════════════════════════
#  求解与延拓
# ═══════════════════════════════════════════════════════════════
def solve_baseflow(Tw_target, verbose=True):
    """
    使用 solve_bvp 求解指定 Tw 的基本流。
    自动从 Tw=1.0 延拓到目标 Tw。
    返回 solve_bvp 的解对象, 或 None (失败时)。
    """
    z = np.linspace(0, ZMAX, N_INITIAL)

    # ── 阶段 1: 先解 Tw=1.0 ──
    if verbose:
        print(f"[1/2] 求解等温基态 Tw=1.0 ...", end=" ", flush=True)

    sol = solve_bvp(vonKarman_ODE, vonKarman_bc(1.0), z,
                    initial_guess(z), tol=TOL, max_nodes=MAX_NODES)

    if not sol.success:
        if verbose:
            print(f"❌\n       {sol.message}")
        return None

    Hinf = sol.sol(ZMAX)[0]
    Fp0  = sol.sol(0)[1]
    Gp0  = sol.sol(0)[3]
    if verbose:
        print(f"✅")
        print(f"       H(∞)={Hinf:.6f}  F'(0)={Fp0:.6f}  G'(0)={Gp0:.6f}")

    if abs(Tw_target - 1.0) < 1e-10:
        return sol

    # ── 阶段 2: 延拓到目标 Tw ──
    is_hot = Tw_target > 1.0
    dTw = 0.002 if is_hot else 0.05
    direction = 1 if is_hot else -1

    if verbose:
        wall_type = "热壁" if is_hot else "冷壁"
        print(f"[2/2] 延拓到 Tw={Tw_target:.4f} (dTw={dTw}, {wall_type}) ...",
              end=" ", flush=True)

    current_guess = sol.sol(z)
    current_tw = 1.0
    result = sol  # 防止空范围导致 UnboundLocalError

    # 构造延拓序列，确保范围非空（处理目标 Tw 接近 1.0 的情况）
    tw_vals = np.arange(1.0 + direction*dTw, Tw_target + direction*dTw*0.5, direction*dTw)
    if len(tw_vals) == 0:
        tw_vals = [Tw_target]

    for tw in tw_vals:
        if abs(tw - current_tw) < 1e-10:
            continue

        result = solve_bvp(vonKarman_ODE, vonKarman_bc(tw), z,
                           current_guess, tol=1e-6, max_nodes=MAX_NODES)

        if result.success:
            current_guess = result.sol(z)
            current_tw = tw
        else:
            # 减半步长
            tw_mid = (current_tw + tw) / 2.0
            result = solve_bvp(vonKarman_ODE, vonKarman_bc(tw_mid), z,
                               current_guess, tol=1e-6, max_nodes=MAX_NODES)
            if result.success:
                current_guess = result.sol(z)
                current_tw = tw_mid
                # 再试原目标
                result = solve_bvp(vonKarman_ODE, vonKarman_bc(tw), z,
                                   current_guess, tol=1e-6, max_nodes=MAX_NODES)
                if result.success:
                    current_guess = result.sol(z)
                    current_tw = tw
                else:
                    if verbose:
                        print(f"❌\n       延拓在 Tw={tw:.4f} 处失败")
                        print(f"       → Boussinesq 近似在此壁温下无自洽解")
                    return None
            else:
                if verbose:
                    print(f"❌\n       延拓在 Tw={tw:.4f} 处失败")
                return None

    if verbose:
        Hinf = result.sol(ZMAX)[0]
        Fp0  = result.sol(0)[1]
        Gp0  = result.sol(0)[3]
        Tp0  = result.sol(0)[5]
        print(f"✅")
        print(f"       H(∞)={Hinf:.6f}  F'(0)={Fp0:.6f}  G'(0)={Gp0:.6f}  T'(0)={Tp0:.6f}")

    return result

# ═══════════════════════════════════════════════════════════════
#  输出 (兼容 Julia readdlm / include)
# ═══════════════════════════════════════════════════════════════
def save_results(sol, Tw):
    """保存数据文件 (3 种格式，供 Julia 选用)"""
    N_out = 2000
    z_out = np.linspace(0, ZMAX, N_out)
    y_out = sol.sol(z_out)

    H_arr  = y_out[0]   # 轴向速度
    F_arr  = y_out[2]   # 径向速度
    G_arr  = y_out[4]   # 周向速度
    T_arr  = y_out[6]   # 温度
    dF_arr = y_out[1]   # F'
    dG_arr = y_out[3]   # G'
    dT_arr = y_out[5]   # T'

    Hinf = float(H_arr[-1])
    Fp0  = float(dF_arr[0])
    Gp0  = float(dG_arr[0])
    Tp0  = float(dT_arr[0])

    tag = f"Tw{Tw:.4f}".replace('.', 'p')

    # ── 格式 1: .dat (DelimitedFiles 读取) ──
    data = np.column_stack([z_out, H_arr, F_arr, G_arr, T_arr, dF_arr, dG_arr, dT_arr])
    fname = f"baseflow_{tag}.dat"
    header = (f"# vonKarman Boussinesq Tw={Tw:.4f} Pr={Pr}\n"
              f"# columns: z H F G T dF dG dT")
    np.savetxt(fname, data, header=header, fmt="%.8e")

    # ── 格式 2: .jl (Julia include 直接赋值) ──
    jl_name = f"baseflow_{tag}.jl"
    with open(jl_name, 'w') as f:
        f.write(f"# Auto-generated by Bone.py\n")
        f.write(f"# Rotating-frame von Karman (Boussinesq)  Tw={Tw:.4f}  Pr={Pr}\n\n")
        # 标量
        f.write(f"const BF_Tw = {Tw:.6f}\n")
        f.write(f"const BF_Pr = {Pr:.6f}\n")
        f.write(f"const BF_Hinf = {Hinf:.10f}\n")
        f.write(f"const BF_Fp0  = {Fp0:.10f}\n")
        f.write(f"const BF_Gp0  = {Gp0:.10f}\n")
        f.write(f"const BF_Tp0  = {Tp0:.10f}\n")
        f.write(f"const BF_N    = {N_out}\n")
        f.write(f"const BF_zmax = {ZMAX:.1f}\n\n")
        # 向量 — 用 readable 格式
        for name, arr in [("BF_z", z_out), ("BF_H", H_arr), ("BF_F", F_arr),
                          ("BF_G", G_arr), ("BF_T", T_arr),
                          ("BF_dF", dF_arr), ("BF_dG", dG_arr), ("BF_dT", dT_arr)]:
            f.write(f"const {name} = [\n")
            for i in range(0, len(arr), 4):
                line_vals = ", ".join(f"{v: .10e}" for v in arr[i:i+4])
                f.write(f"    {line_vals}")
                if i + 4 < len(arr):
                    f.write(",")
                f.write("\n")
            f.write("]\n\n")
        f.write(f'println("✅ Baseflow loaded: Tw={Tw:.4f}, N={N_out}")\n')

    # ── 格式 3: 纯数值 (无注释，readdlm 无需 skipstart) ──
    raw_name = f"baseflow_{tag}_raw.dat"
    np.savetxt(raw_name, data, fmt="%.8e")

    print(f"\n📁 输出文件 (Tw={Tw:.4f}):")
    print(f"   {fname}     — 带注释, Julia: readdlm(..., skipstart=2)")
    print(f"   {raw_name} — 纯数值, Julia: readdlm(...)")
    print(f"   {jl_name}     — Julia include 直接定义变量")

    # ── 剖面图 ──
    if not SAVE_PLOT:
        return

    z_plot = np.linspace(0, 10, 500)
    yp = sol.sol(z_plot)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f'Rotating-frame von Karman Base Flow (Tw={Tw:.4f})', fontsize=13)

    titles = ['F (radial)', 'G (azimuthal)', 'H (axial)', 'T (temperature)']
    ylabels = ['F', 'G', 'H', 'T']
    indices = [2, 4, 0, 6]

    for ax, title, yl, idx in zip(axes.flat, titles, ylabels, indices):
        ax.plot(z_plot, yp[idx], 'b-', lw=1.5)
        ax.set_xlabel('z'); ax.set_ylabel(yl); ax.set_title(title)
        ax.grid(True, alpha=0.3)
        if idx == 6:  # temperature
            ax.axhline(y=1.0, color='gray', ls='--', lw=0.8, label='T(∞)=1')
            ax.axhline(y=Tw, color='red', ls=':', lw=0.8, label=f'Tw={Tw:.3f}')
            ax.legend(fontsize=8)

    plt.tight_layout()
    png_name = f"Bone_profile_Tw{Tw:.2f}.png"
    plt.savefig(png_name, dpi=150)
    plt.close()
    print(f"📊 图片已保存: {png_name}")

# ═══════════════════════════════════════════════════════════════
#  PyCall 接口 — Julia 直接调用
# ═══════════════════════════════════════════════════════════════
def get_baseflow(Tw, verbose=False):
    """
    Julia PyCall 接口:
        using PyCall
        pushfirst!(PyVector(pyimport("sys")."path"), "/path/to/Vonkarmen_bone")
        bone = pyimport("Bone")
        z, H, F, G, T, dF, dG, dT, info = bone.get_baseflow(1.04)

    返回:
        z, H, F, G, T, dF, dG, dT : numpy 数组 (PyCall 自动转为 Julia Vector)
        info : dict, 包含 Hinf, Fp0, Gp0, Tp0, Tw, Pr, success
    """
    sol = solve_baseflow(Tw, verbose=verbose)

    if sol is None:
        return (None,) * 8 + ({"success": False, "Tw": Tw},)

    N_out = 2000
    z_out = np.linspace(0, ZMAX, N_out)
    y_out = sol.sol(z_out)

    Hinf = float(sol.sol(ZMAX)[0])
    Fp0  = float(sol.sol(0)[1])
    Gp0  = float(sol.sol(0)[3])
    Tp0  = float(sol.sol(0)[5])

    info = {
        "success": True,
        "Tw":     float(Tw),
        "Pr":     float(Pr),
        "Hinf":   Hinf,
        "Fp0":    Fp0,
        "Gp0":    Gp0,
        "Tp0":    Tp0,
        "N":      N_out,
        "zmax":   float(ZMAX),
    }

    return (z_out, y_out[0], y_out[2], y_out[4], y_out[6],
            y_out[1], y_out[3], y_out[5], info)

# ═══════════════════════════════════════════════════════════════
#  入口
# ═══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    # ── 批量计算: Tw = 1.00, 1.01, 1.02, 1.03, 1.04 ──
    Tw_list = [1.00, 1.01, 1.02, 1.03, 1.04]

    # 如果命令行指定了 Tw，则只算单个
    if len(sys.argv) > 1:
        Tw_list = [TARGET_TW]

    print("=" * 60)
    print("  动坐标系 von Kármán 基本流求解器 (solve_bvp)")
    print("=" * 60)
    print(f"  Tw 序列: {Tw_list},  Pr = {Pr},  zmax = {ZMAX}")
    print("=" * 60)

    for Tw_val in Tw_list:
        print(f"\n{'─'*60}")
        print(f"  ▶ 开始计算 Tw = {Tw_val:.4f}")
        print(f"{'─'*60}")

        sol = solve_baseflow(Tw_val)

        if sol is not None:
            save_results(sol, Tw_val)
            Hinf = sol.sol(ZMAX)[0]
            Fp0  = sol.sol(0)[1]
            Gp0  = sol.sol(0)[3]
            Tp0  = sol.sol(0)[5]
            print(f"\n  ✅ Tw={Tw_val:.4f} 完成!")
            print(f"     壁面: F'(0)={Fp0:.6f}, G'(0)={Gp0:.6f}, T'(0)={Tp0:.6f}")
            print(f"     远场: H(∞)={Hinf:.6f}")
        else:
            print(f"\n  ❌ Tw={Tw_val:.4f} 求解失败，跳过。")
            if len(Tw_list) == 1:
                sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  全部完成! 共处理 {len(Tw_list)} 个壁温。")
    print(f"{'='*60}")


