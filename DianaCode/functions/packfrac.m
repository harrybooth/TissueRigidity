function phi = packfrac(c,ci,m)

phimin = 0.9;
phimax = 0.999;%99
%ci = 1;
%m = 10;

%phi = phimin + (phimax-phimin)*(c.^m)./(c.^m+ci^m);
phi = phimax - (phimax-phimin)*(c.^m)./(c.^m+ci^m);

end