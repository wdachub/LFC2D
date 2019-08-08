function int = BxInitialFun( y,x1,x2,t,vel,pfun)

if isequal(vel, @velRotation)
    
    if isequal(pfun, @pfunPowerRotationX4)
        int=(1/5).*((-1).*x1+x2).^(-1).*cos(t).^4.*((-1).*(x1+y.*tan(t)).^5+( ...
            x2+y.*tan(t)).^5);
    elseif   isequal(pfun, @pfunPowerRotationX1)
        
        int= (1/2).*(x1+x2).*cos(t)+y.*sin(t);
    end
elseif  isequal(vel, @velIntersection)
    if isequal(pfun, @pfunIntersection)
        int=(1/3).*exp(1).^((-4).*t).*(x1.^2+x1.*x2+x2.^2+3.*y.^2);
    end
    
end

end
