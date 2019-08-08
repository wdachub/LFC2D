function plotComponentsD(addData,hFigure)
% Plots the components of the bounding polygon that are stored in the
% struct 'BoundingPolygonComponents' created by CalcBoundaryD(...) or 
% CalcBoundaryD_Adaptive(...)

Linewidth = 2;

%% Read out required Data
BoundingPolygonComponents = addData.BoundingPolygonComponents;

if isfield(BoundingPolygonComponents,'Surface')
    % For the unadaptive code
    surface = BoundingPolygonComponents.Surface;
else
    surface = BoundingPolygonComponents.RefinedSurface;
end
streakLine1 = BoundingPolygonComponents.StreakLine1;
streakLine2 = BoundingPolygonComponents.StreakLine2;
bwImageSurface = BoundingPolygonComponents.BWimageSurface;

%% Create Plot
if exist('hFigure','var')
    figure(hFigure)
    hold on
else
    figure
    hold on
    axis equal
end
title('Components of \(\partial D\)','interpreter','latex')
grid on
set(gca,'ColorOrderIndex',1)    % Begin with first color of default order
p1 = plot(streakLine1(:,1),streakLine1(:,2),'LineWidth',Linewidth);
p2 = plot(streakLine2(:,1),streakLine2(:,2),'LineWidth',Linewidth);
p3 = plot(surface(:,1),surface(:,2),'LineWidth',Linewidth);
p4 = plot(bwImageSurface(:,1),bwImageSurface(:,2),'LineWidth',Linewidth);
legend([p1 p2 p3 p4],'Streakline1','Streakline2','Surface','Backadvection of Surface','location','best')
xlabel('x')
ylabel('y')
end
