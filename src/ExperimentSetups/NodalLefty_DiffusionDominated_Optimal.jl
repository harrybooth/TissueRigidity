const DN0 = 1.95
const DL0 = 15.
const kN0 = 1.5194850483459551e-6
const kL0 = 7.465310302629284e-7
const kE =  0.00028852347539516876
const kNL = 10.303407667201299
const σN0 = 0.01
const σL0 = 0.00046538852919702296 # alpha_L in Diana?
const Na = 31.6228
const NL = 181.38301582100422
const NE = 12.316816453909029
const mN = 2 # ma in Diana?
const mL = 8 
const mNL = 2
const LN = 19.150715022082355
const s0 = 5. # s0 / 2 as diana scales by dx/2


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
# const α0 = 0.89
const ϕ0 = ϕ(α0)
# const β = 1.
const β = 0.1

const E_star = α_min^(-1) - α0^(-1)

########

const default_mN = mN
const default_mL = mL
const default_mNL = mNL

#########

const λ = sqrt(DN0/kN0)
const c0  = s0/(kN0*λ*(1-exp(-L/λ)))

N_sim = 100000

γ = 0.8
j = 2.

pv_orig = [DN0,DL0,kN0,kL0,kE,kNL,σN0,σL0,Na,NL,NE,LN,s0]

const de_abstol = 1e-10
const de_reltol = 1e-8