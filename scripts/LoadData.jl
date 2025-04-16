data = XLSX.readxlsx(scriptsdirx() * "/SmadCalibrationData.xlsx")["3F_Input"] 

data = XLSX.eachtablerow(data) |> DataFrames.DataFrame

const exp_times = Float64.(data[:,"Time Point"]) * 60;
const exp_t0_wt = exp_times[argmax(data[:,"WT"])]
const exp_t0_cp = exp_times[argmax(data[:,"SLB"])]

const max_exp_wt  = maximum(data[:,"WT"])
const max_exp_cp  = maximum(data[:,"SLB"])

const exp_times_times_norm = exp_times ./ exp_t0_wt;

alpha_data = XLSX.readxlsx(scriptsdirx() * "/SmadCalibrationData.xlsx")["WT_alpha"] 
alpha_data = XLSX.eachtablerow(alpha_data) |> DataFrames.DataFrame

alpha_data_times = [0,30,60,90,120] .* 60

alpha_positions = Float64.(alpha_data[:,"Position"])

const alpha_x = Int.((alpha_positions[1:end]  .+ 15))

const alpha_data_times_norm = alpha_data_times ./ exp_t0_wt

function get_exp_summary_metrics(xmax_data)

    exp_t0_wt
    exp_t0_cp

    wt_xMax = maximum(xmax_data[:,"WT"])
    cp_xMax = maximum(xmax_data[:,"SLB"])

    wt_d0 = xmax_data[end,"WT"] ./ wt_xMax
    cp_d0 = xmax_data[end,"SLB"]  ./ cp_xMax
    lm_d0 = 1.

    xmax_peak_ratio = exp_t0_cp / exp_t0_wt

    xmax_mse = 0.

    alpha_mse = 0.

    return (wt_t0 = exp_t0_wt,cp_t0 = exp_t0_cp,wt_xMax = wt_xMax,cp_xMax = cp_xMax,wt_d0 = wt_d0,cp_d0 = cp_d0,lm_d0 = lm_d0,xmax_peak_ratio = xmax_peak_ratio,xmax_mse = xmax_mse,alpha_mse = alpha_mse)
end