

%% Core workflow:
model_dist = 1; reps = 15;
for mm = 1
    for session = 1
        [hexa_data]     = HX_load_csv([path filenames{mm}], 0, 0);
        [hexa_data_an]  = HX_analyze_session(hexa_data,session);
        if model_dist
            pfm = figure(61); clf;
            for jj=1:reps
                [hexa_model]    = HX_model_session(hexa_data_an,'e-proportional','matchP-shift-spatial',1,0);
                hexa_model_an.sim_reps.ideal(jj,:) = cumsum(sum(hexa_model.ideal,1));
                hexa_model_an.sim_reps.random(jj,:) = cumsum(sum(hexa_model.random,1));
                hexa_model_an.sim_reps.rewards(jj,:) = cumsum(sum(hexa_model.rewards,1));                
            end
            subplot(131); hold off;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.ideal,1),std(hexa_model_an.sim_reps.ideal,1),{'color',[0 0 0]}); hold on;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.random,1),std(hexa_model_an.sim_reps.random,1),{'color',[0.5 0.5 0.5]}); hold on;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.rewards,1),std(hexa_model_an.sim_reps.rewards,1),{'color',[0.5 0 0.16]}); hold on;
            plot(1:size(hexa_model_an.sim_reps.ideal,2),cumsum(sum(hexa_data_an.rewards,1)),'-','color',[1 0 0.33],'linewidth',2); hold on;
            box off;
            ylabel('Cumulative rewards'); xlabel('Unique Port Visits');
            axis([0 size(hexa_model.ideal,2) 0 max(cumsum(sum(hexa_model.ideal,1)))]);
        else
            [hexa_model]    = HX_model_session(hexa_data_an,'e-proportional','matchP-shift-spatial',1,1);
        end
    end
end

%% Filenames
% filenames = {'6PG5_NAc_conc_behav.h5','ML12_NAc_conc_behav.h5','ML13_NAc_conc_behav.h5','ML14_DMS_conc_behav.h5'}
filenames = {'6PG5_NAc_conc_beh.csv'};
path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';

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

