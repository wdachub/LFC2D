function pfun = pfunRotationTimeDep(x, y,t)
theta = atan(t);
x0=x*cos(theta) + y*sin(theta);
y0=-x*sin(theta) + y*cos(theta);
pfun=100*x0.*y0;
end