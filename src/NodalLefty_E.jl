α(E) = 1 / (E + α0^(-1))

ν(N,σ,a,h) = σ*(1/(1+(a/N)^h))

D(D0,ϕ0,ϕ) = (D0/ϕ0)*ϕ

k(k0,ϕ0,β,ϕ) = k0*((β + ϕ0)/(1 - ϕ0))*((1 - ϕ)/(β + ϕ))

σ(σ0,ϕ0,ϕ) = (σ0/(1-ϕ0))*(1-ϕ)

ϕ_tri(α) = 1 - (sqrt(3)/(2*α^2))*(pi/3-2*acos(α)+2*α*sqrt(1-α^2));
ϕ_quad(α) = 1 - (1/(4*α^2))*(pi-4*acos(α)+4*α*sqrt(1-α^2));

ϕ(α) = ϕ_min + ϕ_tri(α)*Int(α>αc_tri) + ϕ_quad(α)*Int(α>αc_quad) 

softplus(x,k) = log(1+exp(k*x)) / k

function nodal_lefty_spatial_diff!(du,u,p,t)

    DN0,DL0,kN0,kL0,kE,kNL,σN0,σL0,Na,NL,NE,mN,mL,mNL,LN,s0 = p

    N = @view u[:,1]
    L = @view u[:,2]
    E = @view u[:,3]
    α = @view u[:,4]

    dN =  @view du[:,1]
    dL =  @view du[:,2]
    dE =  @view du[:,3]
    dα =  @view du[:,4]

    DN = D.(DN0,ϕ0,ϕ.(α))
    DL = D.(DL0,ϕ0,ϕ.(α))
    kN = k.(kN0,ϕ0,β,ϕ.(α))
    kL = k.(kL0,ϕ0,β,ϕ.(α))
    σN = σ.(σN0,ϕ0,ϕ.(α))
    σL = σ.(σL0,ϕ0,ϕ.(α))

    σE = E_star*kE

    h = 1/dx^2

    dN[1] = DN[1]*h*(N[2] - N[1] + (s0/DN[1])) - kN[1]*N[1] - kNL*N[1]*(1/(1+(LN/L[1])^mNL)) + ν(N[1],σN[1],Na,mN)
    dL[1] = DL[1]*h*(L[2] - L[1]) - kL[1]*L[1] + ν(L[1],σL[1],NL,mL)
    dE[1] = ν(N[1],σE,NE,1) - kE*E[1]
    dα[1] =  -(E[1] + α0^(-1))^(-2)*(ν(N[1],σE,NE,1) - kE*E[1])

    @inbounds for j in 2:Nc-1
        dN[j] =  h*(DN[j]*(N[j-1] + N[j+1] - 2*N[j]) + (DN[j+1] - DN[j])*(N[j+1]- N[j])) - kN[j]*N[j]  - kNL*N[j]*(1/(1+(LN/L[j])^2)) + ν(N[j],σN[j],Na,mN)
        dL[j] =  h*(DL[j]*(L[j-1] + L[j+1] - 2*L[j]) + (DL[j+1] - DL[j])*(L[j+1]- L[j])) - kL[j]*L[j] + ν(N[j],σL[j],NL,mL)
        dE[j] = ν(N[j],σE,NE,1) - kE*E[j]
        dα[j] =  -(E[j] + α0^(-1))^(-2)*(ν(N[j],σE,NE,1) - kE*E[j])
    end

    dN[Nc] = DN[Nc]*h*(N[Nc-1] - N[Nc]) - kN[Nc]*N[Nc] - kNL*N[Nc]*(1/(1+(LN/L[Nc])^2)) + ν(N[Nc],σN[Nc],Na,mN)
    dL[Nc] = DL[Nc]*h*(L[Nc-1] - L[Nc]) - kL[Nc]*L[Nc] + ν(L[Nc],σL[Nc],NL,mL)
    dE[Nc] = ν(N[Nc],σE,NE,1) - kE*E[Nc]
    dα[Nc] =  -(E[Nc] + α0^(-1))^(-2)*(ν(N[Nc],σE,NE,1) - kE*E[Nc])

    dα[dα .> 0] .= 0.

end