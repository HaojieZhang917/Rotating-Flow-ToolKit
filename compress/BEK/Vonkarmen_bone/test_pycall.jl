# test_pycall.jl — 验证 Julia → Python → Julia 通路
using PyCall

# 添加 Bone.py 所在目录到 Python path
pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)

# 导入 Bone 模块 (Bone.py)
bone = pyimport("Bone")

# 调用求解器
println("Calling Python solver via PyCall...")
z, H, F, G, T, dF, dG, dT, info = bone.get_baseflow(1.04)

println("\n=== 结果 ===")
println("success = ", info["success"])
println("Tw      = ", info["Tw"])
println("H(∞)    = ", info["Hinf"])
println("F'(0)   = ", info["Fp0"])
println("G'(0)   = ", info["Gp0"])
println("T'(0)   = ", info["Tp0"])
println("N       = ", info["N"])
println()
println("z 类型:  ", typeof(z),  "  长度: ", length(z))
println("H 类型:  ", typeof(H),  "  长度: ", length(H))
println("z[1:3] = ", z[1:3])
println("H[1:3] = ", H[1:3])
