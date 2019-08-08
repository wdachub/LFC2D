function [h,ax] = plotCrossChannelMixer(CrossChannelMixer,xmin,xmax)
% Creates the plot of the Channel defined by the struct 'Channel' in the
% range [xmin xmax].
%
% List of Parameters:
%
% =========================================================================
%   INPUT
%       Name  |         Description
% ------------|------------------------------------------------------------
% CrossChannelMixer 
%             | Struct storing the parameters for cross-channel-mixer. Use
%             | the m-file 'createCrossChannelMixer' for the generation of
%             | this struct.
%     --------|------------------------------------------------------------
%     xmin    | Left x-boundary of the figure
%     --------|------------------------------------------------------------
%     xmax    | Right x-boundary of the figure
% -------------------------------------------------------------------------
%           
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
% ------------|------------------------------------------------------------
%     h       | Handle of the figure
%     --------|------------------------------------------------------------
%     ax      | Handle of the axes 
% -------------------------------------------------------------------------


e = CrossChannelMixer.e;
f = CrossChannelMixer.f;
p = CrossChannelMixer.p;                   
r = CrossChannelMixer.r;      
nCrossChannels = CrossChannelMixer.nCrossChannels;       %number of the cross-channels

upBoundChannel = p + r;       %Upper boundarys of the channels
lowBoundChannel = p - r;      %Lower boundarys of the channels

h = figure('Position',[100 100 1280 350]);
grid on
hold on
%% Plot lines and shading of the walls 
LineWidthWall = 1;
ColorShading = [.83 .83 .83];


fill([xmin,xmin,lowBoundChannel(1),lowBoundChannel(1)], ...
     [f+(f+e)/2, f, f, f+(f+e)/2],...
     ColorShading,...
      'EdgeColor','None')
  
fill([xmin,xmin,lowBoundChannel(1),lowBoundChannel(1)], ...
     [-e-(f+e)/2, -e, -e, -e-(f+e)/2],...
     ColorShading,...
      'EdgeColor','None')

plot([xmin, lowBoundChannel(1)],[f f],'k','LineWidth',LineWidthWall)
plot([xmin, lowBoundChannel(1)],[-e -e],'k','LineWidth',LineWidthWall)
plot([lowBoundChannel(1),lowBoundChannel(1)],[f;f+(f+e)/2],'k','LineWidth',LineWidthWall)
plot([lowBoundChannel(1),lowBoundChannel(1)],[-e;-e-(f+e)/2],'k','LineWidth',LineWidthWall)

for i=1:nCrossChannels-1
     fill([upBoundChannel(i),upBoundChannel(i),lowBoundChannel(i+1),lowBoundChannel(i+1)], ...
          [f+(f+e)/2, f, f, f+(f+e)/2],...
          ColorShading,...
          'EdgeColor','None')
     fill([upBoundChannel(i),upBoundChannel(i),lowBoundChannel(i+1),lowBoundChannel(i+1)], ...
          [-e-(f+e)/2, -e, -e, -e-(f+e)/2],...
          ColorShading,...
          'EdgeColor','None')
    
    plot([upBoundChannel(i),lowBoundChannel(i+1)],[-e,-e],'k','LineWidth',LineWidthWall)
    plot([upBoundChannel(i),lowBoundChannel(i+1)],[f f],'k','LineWidth',LineWidthWall)
    plot([upBoundChannel(i),upBoundChannel(i)],[f;f+(f+e)/2],'k','LineWidth',LineWidthWall)
    plot([upBoundChannel(i),upBoundChannel(i)],[-e;-e-(f+e)/2],'k','LineWidth',LineWidthWall)
    plot([lowBoundChannel(i+1),lowBoundChannel(i+1)],[f;f+(f+e)/2],'k','LineWidth',LineWidthWall)
    plot([lowBoundChannel(i+1),lowBoundChannel(i+1)],[-e;-e-(f+e)/2],'k','LineWidth',LineWidthWall)
    

end

fill([upBoundChannel(end),upBoundChannel(end),xmax,xmax], ...
     [f+(f+e)/2, f, f, f+(f+e)/2],...
     ColorShading,...
     'EdgeColor','None')
fill([upBoundChannel(end),upBoundChannel(end),xmax,xmax], ...
     [-e-(f+e)/2, -e, -e, -e-(f+e)/2],...
      ColorShading,...
      'EdgeColor','None')

plot([upBoundChannel(end),xmax],[f f],'k','LineWidth',LineWidthWall)
plot([upBoundChannel(end),xmax],[-e -e],'k','LineWidth',LineWidthWall)
plot([upBoundChannel(end),upBoundChannel(end)],[f;f+(f+e)/2],'k','LineWidth',LineWidthWall)
plot([upBoundChannel(end),upBoundChannel(end)],[-e;-e-(f+e)/2],'k','LineWidth',LineWidthWall)

axis equal
axis([xmin, ... 
      xmax, ...
      -e-(f+e)/2, ...
      f+(f+e)/2]);

ax = gca;
end