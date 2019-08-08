% Creates a struct 'CrossChannelMixer' that stores all the relevant 
% parameters for the definition of the cross-channel-mixer. 
%
%
%   Sketch of the cross-channel-mixer
%
%             |   |      | |           |     |               |     |
%             | 1 |      |2|           |  3  |               |  n  |
%  -----------|   |------| |-----------|     |------ ... ----|     |-------
%       |
%       f
%       |
%  ............ undisturbed Fluid interface ........ ... ..................   
%       e
%       |
%  -----------|   |------| |-----------|     |------ ... ----|     |-------
%             |   |      | |           |     |               |     |
%             |   |      | |           |     |               |     |
%
%
%
%
% The following fields are contained in the struct 'CrossChannelMixer'
%
% -------------------------------------------------------------------------
%   U       Fluid speed at the interface. Used to create a parabolic
%           velocity-profile with speed zero at the walls and speed U at
%           the interface
% -------------------------------------------------------------------------
%   e       Width of the lower fluid
% -------------------------------------------------------------------------
%   f       Width of the upper fluid
% -------------------------------------------------------------------------
%   eps     Strength of the disturbance. Used to scale all the cross
%           channels by the same magnitude 
% -------------------------------------------------------------------------
%   p       Positions of the centers of the cross channels. For n
%           cross-channels, p must be a 1 x n vector
% -------------------------------------------------------------------------
%   v       Maximum velocity at the center of each cross channel, defined 
%           as scaling of epsilon (v = eps * [scaling1 scaling2 ... ]) 
%           For n cross-channels, v must be a 1 x n
%           vector
% -------------------------------------------------------------------------
%   r       Half the width of the cross channels. For n cross-channels, r
%           must be a 1 x n vector
% -------------------------------------------------------------------------
%   phi     Phase shift between the cross-channels. For n cross-channels,
%           phi must be a 1 x n vector
% -------------------------------------------------------------------------
%   omega   Frequency of the fluid sloshig
% -------------------------------------------------------------------------
%   nCrossChannels
%           Number of cross-channels in this cross-channel-mixer
% -------------------------------------------------------------------------


%% Define parameter
U = 1;
e = 1;
f = 0.7;
epsilon = 1.0;
p = [1 2.3 3 4 5.5]; %[1 2 3 4 5];
v = epsilon * [1 0.5 0.3 0.8 1];
r =  [0.1 0.2 0.1 0.3 0.4];%[0.1   0.1   0.1   0.3   0.1];
phi = pi * [1   2   3   4   7/2];
omega = 4;

%% Check consistency and store number of cross-channels
% Check if the sizes of the vectors are consistent and save the number of
% cross-channels in that case
if isequal(numel(p), numel(v), numel(r), numel(phi))
    nCrossChannels = numel(p);
else
    error('Number of Channels is not consistent!');
end

%% Create Channel
CrossChannelMixer = struct('U',U,...
                           'e',e,...
                           'f',f,...
                           'epsilon',epsilon,...
                           'p', p,...
                           'v', v,...
                           'r', r,...
                           'phi', phi,...
                           'omega',omega,...
                           'nCrossChannels',nCrossChannels);

%% Delete variables for clarity of workspace
clear U e f epsilon p v r phi omega nCrossChannels