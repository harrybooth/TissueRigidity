function s = inhibition_NL(cN,cL,k_NL,ca,m)
%Nodal-Lefty inhibition term

s = k_NL*cN.*(cL.^m)./((cL.^m)+ca^m);

end