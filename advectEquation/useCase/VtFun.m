
function int = VtFun(x1,x2,y1,y2,t,vel,pfun)
%return the exact value of the change rate of cell average at time t. 
%the cell is rectangle [x1,x2]*[y1,y2]
if isequal(vel, @velRotation)

    if isequal(pfun, @pfunPowerRotationX4)
        int=(1/30).*(x1+(-1).*x2).^(-1).*(y1+(-1).*y2).^(-1).*(15.*(x1.^4+(-1) ...
            .*x2.^4).*(y1.^2+(-1).*y2.^2).*cos(t).^4+(-24).*(x1.^5+(-1).* ...
            x2.^5).*(y1+(-1).*y2).*cos(t).^3.*sin(t)+(-45).*(x1.^4+(-1).* ...
            x2.^4).*(y1.^2+(-1).*y2.^2).*cos(t).^2.*sin(t).^2+45.*(x1.^2+(-1) ...
            .*x2.^2).*(y1.^4+(-1).*y2.^4).*cos(t).^2.*sin(t).^2+24.*(x1+(-1).* ...
            x2).*(y1.^5+(-1).*y2.^5).*cos(t).*sin(t).^3+(-15).*(x1.^2+(-1).* ...
            x2.^2).*(y1.^4+(-1).*y2.^4).*sin(t).^4+20.*(x1.^3+(-1).*x2.^3).*( ...
            y1.^3+(-1).*y2.^3).*cos(2.*t).*sin(2.*t));
    elseif isequal(pfun, @pfunPowerRotationX1)
        int=(1/2).*((y1+y2).*cos(t)+(-1).*(x1+x2).*sin(t));
    end
    
elseif isequal(vel, @velIntersection)
    if isequal(pfun, @pfunIntersection)
    int=(-4/3).*exp((-4).*t).*(x1.^2+x1.*x2+x2.^2+y1.^2+y1.*y2+y2.^2);
    end
    
    
end


end



