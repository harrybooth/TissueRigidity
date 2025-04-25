function phi = packfracBernat(x)

alphacrit_tri = 0.8662;
alphacrit_quad = 0.707;

phimin = 0.01818;
%b = 4.1049;
%phimin = 0.9;
%phimax = 1;%0.999
%ci = 1;
%m = 10;
b=1;
phiexpr = b*(1-(sqrt(3)./(2*x.^2)).*(pi/3-2*acos(x)+2*x.*sqrt(1-x.^2)));
%phi = zeros(size(x));

phi_tri = 1-(sqrt(3)./(2*x.^2)).*(pi/3-2*acos(x)+2*x.*sqrt(1-x.^2));
phi_quad = 1-(1./(4*x.^2)).*(pi-4*acos(x)+4*x.*sqrt(1-x.^2));

phi = phimin+0.*(x<alphacrit_quad)+phi_tri.*(x>alphacrit_tri)+phi_quad.*(x>alphacrit_quad);

phiold = phimin+0.*(x<alphacrit_tri)+phiexpr.*(x>alphacrit_tri);

%figure
%plot(phi);hold on;
%plot(phiold)
%hold off
%phi = phimin + (phimax-phimin)*(c.^m)./(c.^m+ci^m);
%phi = phimax - (phimax-phimin)*(c.^m)./(c.^m+ci^m);

end