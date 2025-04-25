function D0 = Dscaling_inv(phi0,D_min,D)
% here we define the reference value D(phi0) at start of simulation, for each particular D-scaling

% for Deff=Dmin+(D0-Dmin)*phi we have:

D0 = D_min+(D-D_min)/phi0;

% 
end