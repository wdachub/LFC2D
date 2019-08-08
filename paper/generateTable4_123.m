clear
addpath '../LFC';
addpath '../useCase';
outputDir = 'output/';
nTests=3;
orders=[2,4,6];
option=[1,1];

curve=[1,0;
    0,0];

disp('A canonical DR in an exponential divergent velocity');
t0=0;
te=1;
nSeg0  = 16;
% option=[1];
vel=@velExpDivergent;
pfun=@pfunExpDivergent; 
for i=1:length(orders)
    order=[2;orders(i);orders(i);orders(i)+10];
%  order=[2,4,4,6];
  expEC(i,:) = compareOrder...
    (curve, vel, pfun,t0,te ,nSeg0,order,option,nTests);
end
format short e
disp(expEC);
writeTestResultsToLatexTable(expEC, nSeg0, orders, 'expDivergent.txt', outputDir);


disp('A non-normal DR with intersecting streaklines');
t0=0;
te = 10;
nSeg0  =64;
vel=@velStreakintersect;
pfun=@pfunStreakintersect; 
for i=1:length(orders)
    order=[2;orders(i);orders(i);orders(i)];
  expEC(i,:) = compareOrder...
    (curve, vel, pfun,t0,te ,nSeg0,order,option,nTests);
end

disp(expEC);
writeTestResultsToLatexTable(expEC, nSeg0, orders, 'nonNormalDR.txt', outputDir);


disp('A normal DR of degree 2');
t0=0;
te= 5/4;
nSeg0  = 64;
vel=@velIntersection;
pfun=@pfunIntersection; 
for i=1:length(orders)
    order=[2;orders(i);orders(i);orders(i)];
  expEC(i,:) = compareOrder...
    (curve, vel, pfun,t0,te ,nSeg0,order,option,nTests);
end
disp(expEC);
writeTestResultsToLatexTable(expEC, nSeg0, orders, 'intersection2.txt', outputDir);

disp('A normal DR of degree 3');
t0=0;
te= 9/4;
nSeg0  =64;
vel=@velIntersection;
pfun=@pfunIntersection; 
for i=1:length(orders)
    order=[2;orders(i);orders(i);orders(i)];
  expEC(i,:) = compareOrder...
    (curve, vel, pfun,t0,te ,nSeg0,order,option,nTests);
end
disp(expEC);
writeTestResultsToLatexTable(expEC, nSeg0, orders, 'intersection3.txt', outputDir);


