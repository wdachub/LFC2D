function plotDonatingRegions(addData,C,varargin)
% Plots the donating regions defined by the data stored in addData which is
% the result of calcFlux2D. Only Regions up to winding numbers +/-3 are 
% plotted 
% The following inputs are optional:
%
% 'RegionOfInterest' | Bounding Polyong to which the flux calculation was
%                    | restricted. Plotted if handed over
% -------------------|-----------------------------------------------------
%          'hFigure' | Plot of the donating regions is added to the figure
%                    | specified by this handle if handle exists

%% Set properties for the plot
%Define colors for the filling of the different regions
color1 = 'r';           %Color for wNo = 1
colorM1 = 'g';          %Color for wNo = -1
color2 = 'c';           %Color for wNo = 2
colorM2 = 'k';          %Color for wNo = -2
color3 = 'y';           %Color for wNo = 3
colorM3 = 'b';          %Color for wNo = -3

%Define the transparency for the filling of the polygons 
facealpha = 0.3;

%% Parse Input
p = inputParser;
addRequired(p,'addData')
addRequired(p,'C')
addParameter(p,'RegionOfInterest',[])
addParameter(p,'hFigure',[])

parse(p,addData,C,varargin{:})

addData = p.Results.addData;
C = p.Results.C;
hFigure = p.Results.hFigure;
region = p.Results.RegionOfInterest;

%% Initialize structures for the creation of the legend
% Matrices store the handles of the filling process and are used to group 
% fillings of the same winding number. Only for the groups legend entries 
% are created
h1= [];
hM1 = [];
h2 = [];
hM2 = [];
h3 = [];
hM3 = [];

% Boolean matrix to track which winding numbers did apear. Used to create
% the corresponding legend entries
showLegend = false(1,6);

%% Create the plot
%Create figure, label the axes and defining the title
if isempty(hFigure)
    figure
    axis equal
else
    figure(hFigure)
end
hold on
grid on
xlabel('x')
ylabel('y')
title('Donating Regions with winding numbers')

%Loop over all stored winding numbers and check in each step what the
%winding number is
for i=1:numel(addData.WindingNumbers)
    
    %Extract current winding number and polygon
    wNo = addData.WindingNumbers(i);
    Poly = cell2mat(addData.DividedPolygons{:,i});
    
    %Check which winding number is present and plot a filled polygon in the
    %corresponding color. Store the handle afterwards in the corresponding
    %matrix and set the corresponding entry in the 'regionEx'-Matrix to
    %true
    if     wNo == 1
        h = fill(Poly(:,1),Poly(:,2),color1);
        set(h,'facealpha',facealpha);
        h1 = [h1 h];
        showLegend(1) = true;
        
    elseif wNo == -1
        h = fill(Poly(:,1),Poly(:,2),colorM1);
        set(h,'facealpha',facealpha)
        hM1 = [hM1 h];
        showLegend(2) = true;
        
    elseif wNo == 2
        h = fill(Poly(:,1),Poly(:,2),color2);
        set(h,'facealpha',facealpha)
        h2 = [h2 h];
        showLegend(3) = true;
        
    elseif wNo == -2
        h = fill(Poly(:,1),Poly(:,2),colorM2);
        set(h,'facealpha',facealpha)
        hM2 = [hM2 h];
        showLegend(4) = true;
        
    elseif wNo == 3
        h = fill(Poly(:,1),Poly(:,2),color3);
        set(h,'facealpha',facealpha)
        h3 = [h3 h];
        showLegend(5) = true;
        
    elseif wNo == -3
        h = fill(Poly(:,1),Poly(:,2),colorM3);
        set(h,'facealpha',facealpha)
        hM3 = [hM3 h];
        showLegend(6) = true;
        
    end
    
end


%% Create the legend
% Check all the handle matrices and if they are not empty, create a group 
% and set the handles as childs. Afterwards the display in the legend must 
% be enabled 
hGroups = [];
if ~isempty(h1)
    gr_h1 = hggroup;
    set(h1,'Parent',gr_h1);
    set(get(get(gr_h1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
    hGroups = [hGroups , gr_h1];
end

if ~isempty(hM1)
    gr_hM1 = hggroup;
    set(hM1,'Parent',gr_hM1);
    set(get(get(gr_hM1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
    hGroups = [hGroups , gr_hM1];
end

if ~isempty(h2)
    gr_h2 = hggroup;
    set(h2,'Parent',gr_h2);
    set(get(get(gr_h2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
    hGroups = [hGroups , gr_h2];
end

if ~isempty(hM2)
    gr_hM2 = hggroup;
    set(hM2,'Parent',gr_hM2);
    set(get(get(gr_hM2,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
    hGroups = [hGroups , gr_hM2];
end

if ~isempty(h3)
    gr_h3 = hggroup;
    set(h3,'Parent',gr_h3);
    set(get(get(gr_h3,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
    hGroups = [hGroups , gr_h3];
end

if ~isempty(hM3)
    gr_hM3 = hggroup;
    set(hM3,'Parent',gr_hM3);
    set(get(get(gr_hM3,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on');
    hGroups = [hGroups , gr_hM3];
end

%Store the legend text for all cases
legendText = {'WindingNumber = 1','WindingNumber = -1',...
              'WindingNumber = 2','WindingNumber = -2',...
              'WindingNumber = 3','WindingNumber = -3'};
          
if ~isempty(region)
    colorRegion = [1 0.73 0];
    h = fill(region(:,1),region(:,2),colorRegion,'DisplayName','Region');
    set(h,'facealpha',facealpha*0.5);
    plot(region(:,1),region(:,2),'Color',0.6*colorRegion)
    
    %Create Legend Entry
    legendText = {legendText{1,:},'Subset Region'};
    showLegend = logical([showLegend,1]);
    hGroups = [hGroups , h];
end          
          
%Display as legend only the parts of 'legendText' for which the
%corresponding regions exist
[~,hIcons] = legend(hGroups,legendText(showLegend));

%Adjust the transparency of the icons
PatchInLegend = findobj(hIcons,'type','patch');
set(PatchInLegend,'facea',facealpha); 

%Plot the surface
plot(C(:,1),C(:,2),'b','Linewidth',2)

end
    
