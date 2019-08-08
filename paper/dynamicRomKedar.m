addpath '../LFC';
addpath '../useCase';
addpath '../other';
addpath '../other/ComplexVisual';
addpath '../other/altmany-export_fig'; 
te=10;
t0=0;
warning('The computation time may longer than a half hour')
order=[2,6,6];
nSeg=80;
err=0.1;
t=linspace(err,2*pi-err,3*nSeg)';
curve=2*[cos(t),sin(t)];
vel  = @velRomKedar;
[DonatingRegion, TimeCurve]= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
% plotGeneratingCurve(DonatingRegion)

plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[2,1,1,1])
export_fig -eps lobe2.eps
% save lobe.mat DonatingRegion TimeCurve
% print(gcf, '-depsc', ['lobe.eps'])
close(gcf);

%%

te=8.0;
t0=0;
p1=[0.75,0]; 
p2=[0.75,1];
curve=[p1;p2];
RKOrder = 6;
order=RKOrder;
nSeg=60;
vel  = @velRomKedar;
[DonatingRegion, TimeCurve]= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
% plotGeneratingCurve(DonatingRegion)

plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,3,1])
export_fig -eps RomKedar1.eps
% save RomKedar1.mat DonatingRegion TimeCurve
% print(gcf,'RomKedar1.eps','-depsc','-r500');
% print(gcf, '-depsc', ['RomKedar1.eps'])
close(gcf);
%%
% zoom in
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion)
axis([0.5,0.9,0.46,1.08])
[~,indM]=min(abs(DonatingRegion.p-1.03));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-1.56));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-1.77));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-1.41));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-1.92));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-2.05));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-2.51));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
[~,indM]=min(abs(DonatingRegion.p-2.9215));
arrow(DonatingRegion.DR(indM,:),DonatingRegion.DR(indM+1,:),'BaseAngle',30);
export_fig -eps RomKedar1L2.eps
% print(gcf, '-depsc', ['RomKedar1L2.eps'])
% print(gcf,'RomKedar1L2.eps','-depsc','-r500');
% close(gcf);
%%

te=8.0;
t0=0;
p1=[1,0]; 
p2=[1,1];
curve=[p1;p2];
order=[2,6,6];
nSeg=60;
vel  = @velRomKedar;
[DonatingRegion, TimeCurve]= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);


plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,2,1])
NumTicks = 6;
set(gca,'YTick',linspace(0,2,NumTicks))
export_fig -eps RomKedar2.eps
% save RomKedar2.mat DonatingRegion TimeCurve
% print(gcf,'RomKedar1L2.eps','-depsc','-r500');
% print(gcf, '-depsc', ['RomKedar2.eps'])