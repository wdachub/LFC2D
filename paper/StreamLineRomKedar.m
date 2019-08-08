clear
addpath '../LFC';
addpath '../useCase';
addpath '../other';
addpath '../other/ComplexVisual';
addpath '../other/altmany-export_fig'; 

t0=0;
te=4;
dt=0.02;
tt=t0:dt:te;
Mov=moviein(length(tt));

for i=1:length(tt)
[xx,yy] = meshgrid(-3:0.02:3,-3:0.02:3);
f = streamFunvelRomKedar(tt(i), xx,yy);
% maxf=max(max(f))-0.3;
% f(f>maxf)=maxf;
contour(xx,yy,f,30,'k');
axis equal tight
Mov(:,i)=getframe;
end
plotDir = 'output/';

fname=[plotDir,'streamlineRomKedar.gif'];
%  movie(M,1,1)
dt1=dt*5;
for i=1:length(Mov)
    [image,map] = frame2im(Mov(i));
    [im,map2]=rgb2ind(image,128);
    if i==1
        imwrite(im,map2,fname,'GIF','WriteMode','overwrite','DelayTime',dt1,'LoopCount',inf);
    else
        imwrite(im,map2,fname,'WriteMode','append','DelayTime',dt1);
    end
end



