function pfun = pfunStrainTimeDep(x,y,t)
%passive function conserved by time dependent strain flow.
pfun=(4*(x*exp(-t^2/2)).^2+4*(y*exp(t^2/2)).^2).^5;
end   