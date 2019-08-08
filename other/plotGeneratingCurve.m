function y =plotGeneratingCurve(DonatingRegion,Narrow)
% plot Donating Region
% nSeg is the number of segments per edge
%
% Narrow is 4D vector or a number, which indicate the number of arrow on
% the LN, streaklineL, preimage of LN, streaklineN.
indL=1;
indN=find(DonatingRegion.p==1);
indPN=find(DonatingRegion.p==2);
indPL=find(DonatingRegion.p==3);
indC=[1,indN,indPN,indPL,length(DonatingRegion.DR)];


linewidth=1;
plot(DonatingRegion.DR(1:indN,1),DonatingRegion.DR(1:indN,2), 'b','Linewidth',linewidth)%LN
hold on
plot(DonatingRegion.DR(indN:indPN,1),DonatingRegion.DR(indN:indPN,2), 'r','Linewidth',linewidth)%streakline N

plot(DonatingRegion.DR(indPL:end,1),DonatingRegion.DR(indPL:end,2), 'r','Linewidth',linewidth)%steakline L
plot(DonatingRegion.DR(indPN:indPL,1),DonatingRegion.DR(indPN:indPL,2), 'k','Linewidth',linewidth)%preimage of LN
hold off
axis equal tight
grid on 

%%
%add arrow
if exist('Narrow','var')==0
    return
end
if length(Narrow)==1
    Narrow=Narrow*ones(1,4);
end
for indNarrow=1:4
    t=linspace(indNarrow-1,indNarrow, Narrow(indNarrow)+2);
    t(end)=[];
    t(1)=[];%delete the endpoint. we don't want to plot arrow on endpoints.
    for indt=1:length(t)
        plotArrow(DonatingRegion,t(indt),indC);
        
    end
end


end
function plotArrow(DonatingRegion,t,indC)

indt=find(ceil([0,1,2,3]-t)==0);

tmp= diff(DonatingRegion.DR,1,1);
edl =sqrt(sum(tmp.^2,2));
cedl=[0;cumsum(edl)];
L=(cedl(indC(indt+1))-cedl(indC(indt)))*(t-indt+1)+cedl(indC(indt));
[~,indM]=min(abs(cedl-L));

if sum(ismember(indC,indM))
    return
else
    arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),20,'BaseAngle',30);
end
end