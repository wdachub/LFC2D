function plotDonatingRegion(DonatingRegion,curve)

%% Divide Polygon
[BoundingPolygon,simpleLoops,connectedComponents,intersections] = dividePolygon(DonatingRegion.DR);

%% Calculate Winding numbers
[wNo,~,~] = polygonCCwindings(BoundingPolygon,connectedComponents);

%% Filter for NaN's in the winding numbers
NaN_wNo = find(isnan(wNo));
wNo(NaN_wNo) = 0;
flag = false;
if ~isempty(NaN_wNo)
    flag = true;
    fprintf('\n%i winding numbers have been changed from NaN to 0\n\n',numel(NaN_wNo))
end



%% Store Results in struct
addData = struct('BoundingPolygon',[real(BoundingPolygon)',imag(BoundingPolygon)'],... %       'BoundingPolygonComponents',BoundingPolygonComponents,...
                 'Intersections',intersections,...
                 'DividedPolygons',cell2table(simpleLoops),...
                 'WindingNumbers',wNo);

plotDonatingRegionsHK(addData,curve)
end