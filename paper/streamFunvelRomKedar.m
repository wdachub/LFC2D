function f = streamFunvelRomKedar(t, x,y)
%return the OVP velocity field

e=0.1;
if exist('T','var')==0

gamma=0.5;
else
    gamma=T;  
    
end

vv=exp(e)./(2.*besseli(0, e));
yv=exp(e.*(cos(t/gamma)-1));
int=integral(@(s)exp(e*cos(s)),0,t/gamma,'AbsTol',1e-14);
xv=1/2.*gamma./yv.*(t/gamma-2.*vv.*exp(-e).*int);


f=-vv*y+1/2*log(((x-xv).^2+(y+yv).^2)./((x-xv).^2+(y-yv).^2))+x.*y.*sin(t/gamma)*e/gamma;
  
end













