addpath '../useCase';

curve=[1,0;
     0,0]; 
vel=@velIntersection;
pfun=@pfunIntersection; 
nSeg=20;
order=[2,4,4,2];
te=1;
t0=0;
option=[0,0];
tic
flux=fluxDR2D(curve, vel, pfun,t0,te ,nSeg,order,option)
eflux=computeFlux(curve(1,:), curve(2,:), t0, te, vel, pfun)
toc
