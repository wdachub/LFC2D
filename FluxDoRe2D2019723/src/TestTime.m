x1 = [0.75;1.2];
x2 = [0.75;-0.2];
C = [0.75*ones(10,1),linspace(-0.2,1.2,10)'];


TolAngle = 15 * pi/180;
TolLength = 5e-2;
T = 2.5;

v = @(t,x) velocityField(t,x);
%%
[BoundingPolygon,BoundingPolygonComponents] = CalcBoundaryD_Adaptive(v,C,T,TolAngle,TolLength);

[Polygon,simpleLoops,connectedComponents,intersections] = dividePolygon(BoundingPolygon);

[wNo,~,~] = polygonCCwindings(BoundingPolygon,connectedComponents);
% Filter for NaN's in the winding numbers
NaN_wNo = find(isnan(wNo));
wNo(NaN_wNo) = 0;
if ~isempty(NaN_wNo)
    fprintf('\n%i Windingnumbers have been changed from NaN to 0\n\n',numel(NaN_wNo))
end

% Store Results in struct
addData = struct('BoundingPolygon',[real(BoundingPolygon)',imag(BoundingPolygon)'],...
                 'DividedPolygons',simpleLoops,...
                 'WindingNumbers',wNo);

plotDonatingRegions(addData,C)
%%
% tic
% [streakline1,streakline2,tau1,tau2] = CalcStreaklinesAdaptive(x1,x2,T,v,TolAngle,TolLength);
% tStreakline = toc
% 
% tic
% [SectionBA,C] = BackadvectSectionAdaptive(C,T,v,TolAngle,TolLength);
% tSectionBA = toc



% figure
% hold on
% plot(streakline1(:,1),streakline1(:,2))
% plot(streakline2(:,1),streakline2(:,2))
% plot(SectionBA(:,1),SectionBA(:,2))