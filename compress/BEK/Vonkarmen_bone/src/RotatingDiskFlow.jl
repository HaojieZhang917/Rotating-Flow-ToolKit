module RotatingDiskFlow

include("NeutralCurveRunner.jl")
include("SutherlandMarching.jl")

using .NeutralCurveRunner: CurveConfig, compute_neutral_curve,
    run_standard_batch, run_parallel_standard_batch,
    validate_curve_file, validate_standard_batch

const LopezBaseflow = NeutralCurveRunner.LopezBaseflow
const LopezStability = NeutralCurveRunner.LopezStability
const NeutralContinuation = NeutralCurveRunner.NeutralContinuation
const CRD_BF = NeutralCurveRunner.CRD_BF

export CurveConfig, compute_neutral_curve
export run_standard_batch, run_parallel_standard_batch
export validate_curve_file, validate_standard_batch
export NeutralCurveRunner, NeutralContinuation
export LopezBaseflow, LopezStability, CRD_BF
export SutherlandMarching

end # module RotatingDiskFlow
