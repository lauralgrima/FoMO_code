

%% Core workflow:

[hexa_data]     = HX_load_hdf5(path , filename, verbose);
[hexa_data_an]  = HX_analyze_session(hexa_data,session);
[hexa_model]    = HX_model_session(hexa_data_an,policy_type,belief_type,plot_out);

%% Filenames
filenames = {'6PG5_NAc_conc_behav.h5','ML12_NAc_conc_behav.h5','ML13_NAc_conc_behav.h5','ML14_DMS_conc_behav.h5'};

path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/hdf5';

%% Belief model types

% 'win-stay' - biased towards staying at current port after reward; visit with no reward explores
% 'matching' - P(rew|port) = sum(rew(port))
% 'match-shift' - P(rew|port) = num_rew +
%           tendency to shift after a success
% 'matchP-shift' - P(rew|port) = num_rew./num_visits +
%           tendency to shift after a success
% 'matchP-shift-local' - P(rew|port) = num_rew(port)./num_visits(port) +
%           tendency to shift after a success
% 'matchP-shift-spatial' - P(rew|port) = num_rew./num_visits +
%           tendency to shift after a success + discounting by distance
% 'matchP-shift-local-spatial' - P(rew|port) = num_rew./num_visits +
%           tendency to shift after a success + discounting by distance
% 'kernel' - P(rew|port) = decaying P(rew) after reward
% 'spatial' - proportional + discount due to distance to port from current location
% 'hazard' - attempt to estimate true hazard P(rew|port,t)
% 'pdf-space' - combined belief about posterior and discounting by distance

%% Policy model types
% 
% 'softmax'
% 'greedy'
% 'e-greedy'
% 'random'
% 'proportional'
% 'e-proportional'

