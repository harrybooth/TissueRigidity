function Dvec = Dscaling(phivec,D_min,D0)
% here effective diffusivity D(phi) is defined

Dvec = D_min+(D0-D_min).*phivec;

end