function s = relay(c,alpha,ca,m)

s = alpha.*(c.^m)./((c.^m)+ca.^m);

end