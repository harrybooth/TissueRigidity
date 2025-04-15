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