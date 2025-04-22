const Nc = 300
const L = 300.
const tissue = range(0,L,length = Nc)
const dx = step(tissue)

const αc_tri = 0.8662;
const αc_quad = 0.707;
const ϕ_min = 0.001

const α_min = 0.7
const α0 = 0.881586152184903
const ϕ0 = ϕ(α0)
const β = 0.1

const E_star = α_min^(-1) - α0^(-1)

########

const default_mN = mN
const default_mL = mL
const default_mNL = mNL

#########

const DN0 = 1.95
const kN0 = 0.5*1e-2
const s0 = 5. # s0 / 2 as diana scales by dx/2?

const λ = sqrt(DN0/kN0)
const c0  = s0/(KN0*λ*(1-exp(-L/λ)))
const N0 = 0.95*c0
const σ_crit = kN0*N0



const σN0 = 1.5*σ_crit
const Na = N0
const Nrelay = (σN0/kN0) + c0

const NL = Nrelay
const kL0 = 1e-4*kN0
const σL0 = 10*kL0 # alpha_L in Diana?
const LN = 2.
const kNL = 1e8 * kN0
const kE =  4*1e-2*kN0
const NE = 0.5*Na


const DL0 = 15.
const mN = 2 # ma in Diana?
const mL = 8 
const mNL = 2


############

N_sim = 100000

γ = 0.8

pv_orig = [DN0,DL0,kN0,kL0,kE,kNL,σN0,σL0,Na,NL,NE,LN,s0]

const de_abstol = 1e-10
const de_reltol = 1e-8