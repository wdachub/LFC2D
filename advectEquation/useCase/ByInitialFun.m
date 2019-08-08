function int = ByInitialFun( x,y1,y2,t ,vel,pfun)
if isequal(vel, @velRotation)
    
    if isequal(pfun, @pfunPowerRotationX4)
        int=(1/5).*((-1).*y1+y2).^(-1).*((-1).*(y1+x.*cot(t)).^5+(y2+x.*cot(t) ...
            ).^5).*sin(t).^4;
    elseif   isequal(pfun, @pfunPowerRotationX1)
        int=x.*cos(t)+(1/2).*(y1+y2).*sin(t);
    end
elseif  isequal(vel, @velIntersection)
    if isequal(pfun, @pfunIntersection)
        int=(1/3).*exp(1).^((-4).*t).*(3.*x.^2+y1.^2+y1.*y2+y2.^2);
    end
end
end
