function int = ByVel( x,y1,y2,t,vel)
%return the face integral of x component of velocity with respect to  y. 

if isequal(vel, @velRotation)
int=1/2*(-y1 - y2);
elseif isequal(vel, @velIntersection)
int=pi*(y1 + y2) + x;
end


end
