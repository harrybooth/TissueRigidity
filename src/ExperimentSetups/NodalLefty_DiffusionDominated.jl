const DN0 = 1.95
const DL0 = 15.
const kN0 = 1.25*1e-6
const kL0 = 7.5*1e-7
const kE =  5*1e-4
const kNL = 10.
const σN0 = 1e-2
const σL0 = 12*1e-3 # alpha_L in Diana?
const Na = 31.6228
const NL = 200.
const NE = 0.75*Na
const mN = 2 # ma in Diana?
const mL = 8 
const mNL = 2
const LN = 15.8114
const s0 = 5. # s0 / 2 as diana scales by dx/2?

############

const Nc = 300
const L = 300.
const tissue = range(0,L,length = Nc)
const dx = step(tissue)

const αc_tri = 0.8662;
const αc_quad = 0.707;
const ϕ_min = 0.001;

const α_min = 0.7
const α0 = 0.881586152184903
const ϕ0 = ϕ(α0)
# const β = 1.
const β = 0.1

const E_star = α_min^(-1) - α0^(-1)

########

const default_mN = mN
const default_mL = mL
const default_mNL = mNL

#########

N_sim = 10.

γ = 0.8
j = 2.

pv_orig = [DN0,DL0,kN0,kL0,kE,kNL,σN0,σL0,Na,NL,NE,LN,s0]

const de_abstol = 1e-10
const de_reltol = 1e-8