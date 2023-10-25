

%% Core workflow:
model_dist = 1; reps = 20; clear hexa_model_an;
for mm = 1
    for session = 1
        [hexa_data]     = HX_load_csv([path filenames{mm}], 0, 1);
        [hexa_data_an]  = HX_analyze_session(hexa_data,session,1);
        if model_dist
            pfm = figure(61); clf;
            for jj=1:reps
                [hexa_model]    = HX_model_session(hexa_data_an,'e-proportional','match-shift-spatial',1,0);
                hexa_model_an.sim_reps.ideal(jj,:)      = cumsum(sum(hexa_model.ideal,1));
                hexa_model_an.sim_reps.random(jj,:)     = cumsum(sum(hexa_model.random,1));
                hexa_model_an.sim_reps.rewards(jj,:)    = cumsum(sum(hexa_model.rewards,1));                
                hexa_model_an.sim_reps.Sideal(jj,:)     = hexa_model.slope.ideal;
                hexa_model_an.sim_reps.Srandom(jj,:)    = hexa_model.slope.random;
                hexa_model_an.sim_reps.Srewards(jj,:)   = hexa_model.slope.model;    
                hexa_model_an.sim_reps.visits(jj,:)     = sum(hexa_model.visits,2);
            end
            subplot(131); hold off;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.ideal,1),2.*std(hexa_model_an.sim_reps.ideal,1),{'color',[0 0 0]}); hold on;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.random,1),2.*std(hexa_model_an.sim_reps.random,1),{'color',[0.5 0.5 0.5]}); hold on;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.rewards,1),2.*std(hexa_model_an.sim_reps.rewards,1),{'color',[0.5 0 0.16]}); hold on;
            plot(1:size(hexa_model_an.sim_reps.ideal,2),cumsum(sum(hexa_data_an.rewards,1)),'-','color',[1 0 0.33],'linewidth',2); hold on;
            box off;
            ylabel('Cumulative rewards'); xlabel('Unique Port Visits');
            axis([0 size(hexa_model.ideal,2) 0 max(cumsum(sum(hexa_model.ideal,1)))]);

            subplot(132); hold off;
            shadedErrorBar(hexa_model.slope.x,mean(hexa_model_an.sim_reps.Sideal,1),2.*std(hexa_model_an.sim_reps.Sideal,1),{'color',[0 0 0]}); hold on;
            shadedErrorBar(hexa_model.slope.x,mean(hexa_model_an.sim_reps.Srandom,1),2.*std(hexa_model_an.sim_reps.Srandom,1),{'color',[0.5 0.5 0.5]}); hold on;
            shadedErrorBar(hexa_model.slope.x,mean(hexa_model_an.sim_reps.Srewards,1),2.*std(hexa_model_an.sim_reps.Srewards,1),{'color',[0.5 0 0.16]}); hold on;
            plot(hexa_model.slope.x,hexa_model.slope.mouse,'-','color',[1 0 0.33],'linewidth',2); hold on;
            box off;
            ylabel('Local Income'); xlabel('Session Time');
            axis([0 size(hexa_model.ideal,2) 0.015 0.075]);

            subplot(133); hold off;
            shadedErrorBar(1:6,mean(hexa_model_an.sim_reps.visits,1),2.*std(hexa_model_an.sim_reps.visits,1),{'color',[0.5 0 0.16]}); hold on;
            plot(1:6,sum(hexa_data_an.visits,2),'o-','color',[1 0 0.33],'linewidth',2);
            xlabel('Total visits'); ylabel('Port');
            axis([0 7 0 max(sum(hexa_data_an.visits,2))*1.5]); 
            box off;

        else
            [hexa_model]    = HX_model_session(hexa_data_an,'e-proportional','match-shift',1,1);
        end
    end
end

%% Filenames
% filenames = {'6PG5_NAc_conc_behav.h5','ML12_NAc_conc_behav.h5','ML13_NAc_conc_behav.h5','ML14_DMS_conc_behav.h5'}
filenames = {'6PG5_NAc_conc_beh.csv'};
path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';

%% Calculate a transition matrix from visits matrix
tmp = find(sum(hexa_model.visits,1)==1);
[~,visit_list] = max(hexa_model.visits(:,tmp),[],1);

[trans_mat] = HX_ComputeTransitionMatrix(visit_list,25);


tmp = find(sum(hexa_data_an.visits,1)==1);
[~,visit_list_data] = max(hexa_data_an.visits(:,tmp),[],1);

[trans_mat_data] = HX_ComputeTransitionMatrix(visit_list_data,26);


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

