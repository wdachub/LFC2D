addpath '../test';
addpath '../LFC';
addpath '../useCase';
addpath '../other';
addpath '../other/ComplexVisual';
addpath '../other/altmany-export_fig';

%initial data

order=[2,6,6];
plotDir = 'output/';


%%
% % ExpDivergent
t0=0;
te=1.0;
nSeg=8;
p1=[0,0];
p2=[1,0];
curve=[p2;p1];
vel = @velExpDivergent;
DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,1,2])
NumTicks = 6;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTick',linspace(-1,0,NumTicks))
export_fig -eps expDivergent.eps
% print(gcf, '-depsc', [plotDir,'expDivergent.eps'])
close(gcf);

%%
% % DR with degree 2
t0=0;
te=1.25;
nSeg=16;
vel  = @velIntersection;
p1=[0,0];
p2=[1,0];
curve=[p1;p2];
DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,3,0,0])
export_fig -eps intersection2.eps
% print(gcf, '-depsc', [plotDir,'intersection2.eps'])
close(gcf);

%%
%DR with degree 3
t0=0;
te=2.25;
nSeg=30;
vel  = @velIntersection;
p1=[0,0];
p2=[1,0];
curve=[p1;p2];
DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,3,0,0])
export_fig -eps intersection3.eps
% print(gcf, '-depsc',[plotDir,'intersection3.eps'])
close(gcf);


%%
%intersect_streakline with pathline from given points
te= 10;
t0=0;
nSeg=16;
p1=[0,0];
p2=[1,0];
curve=[p1;p2];
vel = @velStreakintersect;
DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,1,2])
hold on
indL=1;
indN=find(DonatingRegion.p==1);
indPN=find(DonatingRegion.p==2);
indPL=find(DonatingRegion.p==3);

plot(p1(1),p1(2),'bs','MarkerSize', 8,'Linewidth',1.5);%L
plot(p2(1),p2(2),'bs','MarkerSize', 8,'Linewidth',1.5);%N
% plot([p1(1),p2(1)],[p1(2),p2(2)],'b','Linewidth',2);
plot(DonatingRegion.DR(indPL,1),DonatingRegion.DR(indPL,2),'kd','MarkerSize', 8,'Linewidth',1.5);%preimage of L
plot(DonatingRegion.DR(indPN,1),DonatingRegion.DR(indPN,2),'kd','MarkerSize', 8,'Linewidth',1.5);%preimage of N
% plot(DonatingRegion.DR([indPL,indPN],1),DonatingRegion.DR([indPL,indPN],2),'k','Linewidth',2);
export_fig -eps streaklineDR1.eps
% print(gcf, '-depsc', [plotDir,'streaklineDR1.eps'])
close(gcf);

%%
%algebraic decomposition
te= 10;
t0=0;
k=5;
nSeg=16;
p1=[0,0];
p2=[1,0];
curve=[p1;p2];
vel = @velStreakintersect;
DonatingRegion= donatingRegion(curve, vel,t0,k ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,1,2])
hold on
indL=1;
indN=find(DonatingRegion.p==1);
indPN=find(DonatingRegion.p==2);
indPL=find(DonatingRegion.p==3);

plot(p1(1),p1(2),'bs','MarkerSize', 8,'Linewidth',1.5);%L
plot(p2(1),p2(2),'bs','MarkerSize', 8,'Linewidth',1.5);%N
plot(DonatingRegion.DR(indPL,1),DonatingRegion.DR(indPL,2),'kd','MarkerSize', 8,'Linewidth',1.5);%preimage of L
plot(DonatingRegion.DR(indPN,1),DonatingRegion.DR(indPN,2),'kd','MarkerSize', 8,'Linewidth',1.5);%preimage of N
export_fig -eps streaklineDRD1.eps
% print(gcf,'streaklineDRD1.eps','-depsc','-r800');
% print(gcf, '-depsc', [plotDir,'streaklineDRD1.eps'])
close(gcf);
% axis([-8 8 -8 5]);
%%
nSeg=16;
curve=[p1;p2];
DonatingRegion= donatingRegion(curve, vel,te-k,te ,nSeg,order,[1,1]);

dt = (te-t0)/nSeg;
%compute the DR.
DonatingRegion.DR=flowmap(DonatingRegion.DR, k,t0 , vel, -dt);

plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,1,2])
hold on
indL=1;
indN=find(DonatingRegion.p==1);
indPN=find(DonatingRegion.p==2);
indPL=find(DonatingRegion.p==3);

plot(DonatingRegion.DR(indL,1),DonatingRegion.DR(indL,2),'bs','MarkerSize', 8,'Linewidth',1.5);% L under inverse flow map 
plot(DonatingRegion.DR(indN,1),DonatingRegion.DR(indN,2),'bs','MarkerSize', 8,'Linewidth',1.5);% N under inverse flow map
plot(DonatingRegion.DR(indPL,1),DonatingRegion.DR(indPL,2),'kd','MarkerSize', 8,'Linewidth',1.5);%preimage of L
plot(DonatingRegion.DR(indPN,1),DonatingRegion.DR(indPN,2),'kd','MarkerSize', 8,'Linewidth',1.5);%preimage of N
export_fig -eps streaklineDRD2.eps
% print(gcf,'streaklineDRD2.eps','-depsc','-r500');
% print(gcf, '-depsc', [plotDir,'streaklineDRD2.eps'])
close(gcf);

%%

te= 10;
t0=0;
tt=t0:dt:te;
x=0:dt:1;
nseg=100;

DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
indL=1;
indN=find(DonatingRegion.p==1);
indPN=find(DonatingRegion.p==2);
indPL=find(DonatingRegion.p==3);
indC=[1,indN,indPN,indPL,length(DonatingRegion.DR)];
plot(DonatingRegion.DR(indPL:end,1),DonatingRegion.DR(indPL:end,2), 'r','Linewidth',1.5)%steakline L
hold on
% SlefxSelf(pts)
% x0=-6.257410859208838;%pathline from streaklines intersect point
% y0=0.258932595953038;
x0=-0.5888;%pathline from self-intersection of a streakline 
y0=-2.3082;
% another intersection
%x0=-2.126023834681789;
%y0=-3.897341595662277;
plot(x0,y0,'bo','Linewidth',2)
plot(0,0,'bs','Linewidth',2)
tt=linspace(te,t0,nseg);
path(1,:)=-(y0+1)*sin(tt)+(5+x0)*cos(tt)+tt-5;
path(2,:)=-(y0+1)*-cos(tt)+(5+x0)*sin(tt)-1;
plot(path(1,:),path(2,:),'k');
axis equal tight
%print(gcf, '-depsc', [plotDir,'DR+pathline.eps'])
export_fig -eps DR+pathline.eps
close(gcf);
%%
plot(DonatingRegion.DR(indN:indPN,1),DonatingRegion.DR(indN:indPN,2), 'r','Linewidth',1.5)%steakline L
hold on
plot(DonatingRegion.DR(indPL:end,1),DonatingRegion.DR(indPL:end,2), 'r','Linewidth',1.5)%steakline L
x0=-6.257410859208838;%pathline from streaklines intersect point
y0=0.258932595953038;
plot(x0,y0,'bo','Linewidth',2)
tt=linspace(7,t0,nseg);
path(1,:)=-(y0+1)*sin(tt)+(5+x0)*cos(tt)+tt-5;
path(2,:)=-(y0+1)*-cos(tt)+(5+x0)*sin(tt)-1;
plot(path(1,:),path(2,:),'k');
plot(0,0,'bs','Linewidth',2)
plot(1,0,'bs','Linewidth',2)
axis equal tight
hold off
%print(gcf, '-depsc', [plotDir,'DR+pathline 2.eps'])
export_fig -eps DR+pathline2.eps
%close(gcf);

% x0=-1.8;%pathline from region with degree=0
% y0=-3.897341595662277;
% plot(x0,y0,'ro')
% tt=linspace(te,t0,nseg);
% path(1,:)=-(y0+1)*sin(tt)+(5+x0)*cos(tt)+tt-5;
% path(2,:)=-(y0+1)*-cos(tt)+(5+x0)*sin(tt)-1;
% plot(path(1,:),path(2,:),'--r');

%%
%RomKedar
te=8;
t0=0;
p1=[0.75,-1]; 
p2=[0.75,1];
curve=[p1;p2];
order=[2,6,6,2];
nSeg=30;
vel  = @velRomKedar;

DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,1,2])
export_fig -eps RomKedar1.eps
print(gcf, '-depsc', [plotDir,'RomKedar1.eps'])
save RomKedar1.mat DonatingRegion

close(gcf);



%%
te=8;
t0=0;
p1=[1,-1]; 
p2=[1,1];
curve=[p1;p2];
order=[2,6,6,2];
nSeg=40;
vel  = @velRomKedar;

DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
plotDonatingRegion(DonatingRegion,curve)
plotGeneratingCurve(DonatingRegion,[1,2,1,2])
save RomKedar2.mat DonatingRegion
export_fig -eps RomKedar2.eps
print(gcf, '-depsc', [plotDir,'RomKedar2.eps'])
close(gcf);



%%
%movie
%the evolution of DR with given velocity field
t0=0;
te=10;
dt=0.05;
tt=t0:dt:10;
p1=[0,0];
p2=[1,0];
curve=[p1;p2];
Mov=moviein(length(tt));
vel  = @velStreakintersect;
nSeg=8;
DonatingRegion= donatingRegion(curve, vel,t0,te ,nSeg,order,[1,1]);
% plotGeneratingCurve(DonatingRegion)
% 
% axis([-8 9 -8 5]);

RKOrder=4;
Mov(:,1)=getframe;
indL=1;
indN=find(DonatingRegion.p==1);
indPN=find(DonatingRegion.p==2);
indPL=find(DonatingRegion.p==3);
% ln=pts(:,Vertices(3):end);
ln=[-0.05,0;1.05,0];%because there is numerical error. streakline may not intersect
% 
 st1=DonatingRegion.DR(indN:indPN,:);
timel=DonatingRegion.DR(indPN:indPL,:);
st2=DonatingRegion.DR(indPL:end,:);

for i=1:length(tt)-1
    % delete the points which have past LN.
    st1 = RungeKutta(st1, tt(i), tt(i+1), vel, RKOrder);
    st2 = RungeKutta(st2, tt(i), tt(i+1), vel, RKOrder);
    timel = RungeKutta(timel, tt(i), tt(i+1), vel, RKOrder);
    [xi,yi,ii]=polyxpoly(st1(:,1), st1(:,2), ln(:,1), ln(:,2));
    if ~isempty(xi)
        st1=[[xi(1),yi(1)];st1(ii(1,1)+1:end,:)];
    end
    [xi,yi,ii]=polyxpoly(st2(:,1), st2(:,2), ln(:,1), ln(:,2));
    if ~isempty(xi)
        st2=[st2(1:ii(end,1),:);[xi(end),yi(end)]];
    end
    plot([1,0],[0,0],'-');
    hold on
    
    if i~=length(tt)-1
        plot(st1(:,1),st1(:,2),'-*');
        plot(st2(:,1),st2(:,2),'-*');
    end
    plot(timel(:,1),timel(:,2),'-d');
    axis([-8 8 -8 5]);
    hold off
    
    Mov(:,i+1)=getframe;
end


fname=[plotDir,'intersect_streaklne.gif'];
%  movie(M,1,1)
dt1=dt/2;
for i=1:length(Mov)
    [image,map] = frame2im(Mov(i));
    [im,map2]=rgb2ind(image,128);
    if i==1
        imwrite(im,map2,fname,'GIF','WriteMode','overwrite','DelayTime',dt1,'LoopCount',inf);
    else
        imwrite(im,map2,fname,'WriteMode','append','DelayTime',dt1);
    end
end

