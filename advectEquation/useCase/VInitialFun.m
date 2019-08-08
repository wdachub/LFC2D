function int = VInitialFun( x1,x2,y1,y2,t,vel,pfun)
%return the exact value of cell average at time t.
%the cell is rectangle [x1,x2]*[y1,y2]


if isequal(vel, @velRotation)
    
    if isequal(pfun, @pfunPowerRotationX4)
        int=(1/30).*(x1+(-1).*x2).^(-1).*(y1+(-1).*y2).^(-1).*(6.*(x1.^5+(-1) ...
            .*x2.^5).*(y1+(-1).*y2).*cos(t).^4+15.*(x1.^4+(-1).*x2.^4).*( ...
            y1.^2+(-1).*y2.^2).*cos(t).^3.*sin(t)+15.*(x1.^2+(-1).*x2.^2).*( ...
            y1.^4+(-1).*y2.^4).*cos(t).*sin(t).^3+6.*(x1+(-1).*x2).*(y1.^5+( ...
            -1).*y2.^5).*sin(t).^4+5.*(x1.^3+(-1).*x2.^3).*(y1.^3+(-1).*y2.^3) ...
            .*sin(2.*t).^2);
    elseif   isequal(pfun, @pfunPowerRotationX1)
        
        int=(1/2).*((x1+x2).*cos(t)+(y1+y2).*sin(t));
    end
elseif  isequal(vel, @velIntersection)
    if isequal(pfun, @pfunIntersection)
        int=(1/3).*exp(1).^((-4).*t).*(x1.^2+x1.*x2+x2.^2+y1.^2+y1.*y2+y2.^2);
    end
    
    
end

end

