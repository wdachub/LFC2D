function [errConv, errRich] =compareOrder...
    (curve, vel, pfun,t0,te ,nSeg0,order,option,nTests)

r    = 2; 
fluxSum = zeros(1, nTests);
errConv = zeros(1, 2*nTests-1);
errRich = zeros(1, 2*nTests-3);

exactFlux=0;
for k=1:size(curve,1)-1
exactFlux=exactFlux+computeFlux(curve(k,:), curve(k+1,:), t0, te, vel, pfun);
end

for k=1:nTests
    nSeg = nSeg0*r^(k-1);

    fluxSum(:,k) = fluxDR2D(curve, vel, pfun,t0,te ,nSeg(1),order,option);
    
    err = abs((fluxSum(:,k)-exactFlux)/exactFlux);
%     err= abs((fluxSum(:,k)-exactFlux));
    errConv(2*k-1) = norm(err,inf);
    % Richardson error estimate
    if k>1
        errRich(2*k-3) = max(abs(fluxSum(:,k)-fluxSum(:,k-1)));
    end
    if k>2
        errRich(2*k-4) = log(errRich(2*k-5)/errRich(2*k-3)) / log(r);
    end
end
cid = 2:2:2*nTests-2;
errConv(cid) = log(errConv(cid-1)./errConv(cid+1))/log(r);

