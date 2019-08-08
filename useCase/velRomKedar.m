function vel = velRomKedar(t, xy, T)

x = xy(:,1);
y = xy(:,2);


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


rm=(x-xv).^2+(y-yv).^2;
rp=(x-xv).^2+(y+yv).^2;
vel = [-((y-yv)./rm-(y+yv)./rp)-vv+e/gamma*x.*sin(t/gamma);
      (x-xv).*(1./rm-1./rp)-e/gamma*y.*sin(t/gamma)];

end













