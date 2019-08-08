function uCrossChannelMixer = VelocityFieldCrossChannelMixer(t,Pos,CrossChannelMixer)
% Definition of the velocity field in the cross-channel-mixer defined by 
% the parameters of the struct 'CrossChannelMixer' created by the m-file 
% 'createCrossChannelMixer'.
% This definition of the velocity field is suitable for the use as 'u' in 
% 'fluxLagrangeSteadySurface2D'
%
% List of Parameters:
%
% =========================================================================
%   INPUT
%       Name  |         Description
% ------------|------------------------------------------------------------
%      t      | Scalar value representing the current time
%   ----------|------------------------------------------------------------
%      Pos    | Positions for that the velocity shall be evaluated. Must 
%             | have the following form:
%             |          [x1
%             |           x2
%             |           ...  
%             |     Pos = xn
%             |           y1
%             |           y2
%             |           ...
%             |           yn]
%   ----------|------------------------------------------------------------
% CrossChannelMixer 
%             | Struct storing the parameters for cross-channel-mixer. Use
%             | the m-file 'createCrossChannelMixer' for the generation of
%             | this struct.
% -------------------------------------------------------------------------
%
%
% =========================================================================
%   OUTPUT
%       Name  |         Description
% ------------|------------------------------------------------------------
% uCrossChannelMixer
%             | The resulting velocities for the given time at the given
%             | positions, stored as follows
%             |                        [vx(t,x1,y1)
%             |                         vx(t,x2,y2)
%             |                         ...  
%             |    uCrossChannelMixer = vx(t,xn,yn)
%             |                         vy(t,x1,y1)
%             |                         vy(t,x2,y2)
%             |                         ...
%             |                         vy(t,xn,yn)]
% -------------------------------------------------------------------------


half = numel(Pos)/2;

%% Get parameter for CrossChannelMixer
U = CrossChannelMixer.U;
e = CrossChannelMixer.e;
f = CrossChannelMixer.f;
p = CrossChannelMixer.p;
v = CrossChannelMixer.v;
r = CrossChannelMixer.r;  
phi = CrossChannelMixer.phi;
omega = CrossChannelMixer.omega;
CrossnChannels = CrossChannelMixer.nCrossChannels;      %number of the cross-channels

%% Steady velocity Field of the channel
uChannel = zeros(size(Pos));
uChannel(1:half) =        -U/e/f*Pos(half+1:end).^2 ...
            + U*(1/e-1/f)*Pos(half+1:end)    ...
            + U;
        
% Set negative values to 0 (outside the channel)     !!!CHECK IF REQUIRED!!! 
uChannel(uChannel<0)=0;
        
%% Unsteady velocity Field
upBoundChannel = p + r;       %Upper boundarys of the channels
lowBoundChannel = p - r;      %Lower boundarys of the channels

% Step 1: Determine for each point of Pos which lower boundaries are 
% smaller than the x component and which upper boundaries are bigger. If 
% for a point the lower boundary of the j-th cross-channel is smaller AND 
% the upper boundary of the same channel is bigger, then the point lies in 
% the impact area of this cross-channel and is thus affected by the flow of
% this cross-channel.  
% Create the logic matrix 'LogicPointChannel' such that the entry i,j = 1 
% if the i-th point of Pos is affected by the j-th channel
smaller = repmat(Pos(1:half),1,CrossnChannels)-repmat(lowBoundChannel,half,1) > 0;
bigger  = repmat(Pos(1:half),1,CrossnChannels)- repmat(upBoundChannel,half,1) < 0;

LogicPointChannel = smaller .* bigger;


% Step 2: calculate for each point the velocity due to the cross-channels
A = v./(r.^2).*cos(omega*t+phi);    %Vector of values that are constant for each cross-channel

uCrossChannels = ((Pos(1:half) - LogicPointChannel*p').^2 -(LogicPointChannel*r').^2) .* (LogicPointChannel*A');

%% Resulting velocity field
uCrossChannelMixer = uChannel + [zeros(half,1);uCrossChannels];
end