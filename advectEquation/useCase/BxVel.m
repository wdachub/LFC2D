function int = BxVel( y,x1,x2,t,vel)
%return the face integral of y component of the velocity with respect to x.

if isequal(vel, @velRotation)
int=(x1+x2)/2;
elseif isequal(vel, @velIntersection)
int=-pi*(x1 + x2) + y;
end

end