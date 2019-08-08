clear
addpath '../LFC';
addpath '../useCase';
addpath '../other';
addpath '../other/ComplexVisual';
addpath '../other/altmany-export_fig';

outputDir = 'output/';
warning('The computation time may longer than 2 hours')
te=8.0;
t0=0;
p1=[0.75,0]; 
p2=[0.75,1];
curve=[p1;p2];
order=[2,4,4,2];
vel  = @velRomKedar;
pfun=@pfunConstant;
option=[1,1];
nSeg=20*2.^(0:4);
for indnSeg=1:length(nSeg)
    tic
flux(indnSeg) = fluxDR2D...
    (curve, vel, pfun,t0,te ,nSeg(indnSeg),order,option)
toc
end

save rkfluxorder4.mat flux
%%

for i=1:length(flux)-1
   reslut(2*i-1)=abs(flux(i+1)-flux(i));
end

for i=1:length(flux)-2
   reslut(2*i)=log(abs(reslut(2*i-1))/abs(reslut(2*i+1)))/log(2);
end
reslut
orders=[6];
nSeg0=10;
writeTestResultsToLatexTable(reslut, nSeg0, orders, 'dynamicRKA.txt',outputDir);