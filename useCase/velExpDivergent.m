function vel = velExpDivergent(time, xy, T)
% return the velocity of a divergent velocity of exponential functions.
% This velocity and the scalar $f=exp(-(x+y+t))$
%  satisfy the scalar conservation law.

x = xy(:,1);
y = xy(:,2);
vel = [exp(x)-1-time, exp(y)+time];

