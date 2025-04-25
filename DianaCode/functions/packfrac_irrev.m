function phi_tplusdt = packfrac_irrev(phi_t,c,ci,m)
%updates packing fraction in an irreversible manner

%phi_new = packfrac(c,ci,m);
phi_new = 1-packfracBernat(c);

phi_tplusdt = phi_t;
phi_tplusdt(phi_new>phi_t) = phi_new(phi_new>phi_t);

% figure(100)
% plot(phi_tplusdt-phi_t);
% axis square

end