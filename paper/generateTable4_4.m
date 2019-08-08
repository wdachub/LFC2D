% This code computes the CPU time and relative errors of flux by LFC method.
% The detailed setting is described in table 4.4. 
% One can compare the CPU time produced by this file with the CPU time
% produced by the file
%  ../FluxDoRe2D2019723/examples/generateTable4 4.m

clear
addpath '../LFC';
addpath '../useCase';
outputDir = 'output/';

nseg0 = 10;
nseg=nseg0.*2.^(1:4);
p1=[0,0];
p2=[1,0];
curve = [p1;p2];
t0=0;
vel=@velIntersection;%veocity field (x+2*pi*y, -2*pi*x+y)
pfun=@pfunIntersection;%scalar field (x^2+y^2)*e^(-4*t)
order=6;
te=2.25;
exactflux=computeFlux(p1, p2, t0, te, vel, pfun);

for cn=1:length(nseg)
    tic
    option=[1,1];
    numflux=fluxDR2D(curve, vel, pfun,t0,te ,nseg(cn),order,option);
    ROutPut(cn)=abs(numflux-exactflux)/abs(exactflux);%relative error
    toc

end
format short e
ROutPut
 formSt = '%2.2e & %2.2e& %2.2e & %2.2e& %2.2e & ';

fid = fopen('table4_4.txt', 'w');
 fprintf(fid, formSt,ROutPut);
 fclose(fid);
