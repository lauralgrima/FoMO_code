function [trans_r2, income_r2, visits_for_LL, rewards_for_LL, trans_r2_rand] = HX_model_session_OIO(visit_matrix,cost_per_port,rew_sched,income,prior)
% Creating a simplified version of model code to allow optimization of
% alpha as a function of tau1 and tau2

    epsilon = 0.05;

    portq = sum(rew_sched,2);

    % sampling rate is now set to 1 Hz
    frame_rate = 1;
    
    % for testing
    % alpha_params = [0.05 0 0 1 1 0.1];

    % Make the BigBoi reward environment
    max_reward = sum(sum(rew_sched));

    max_tsteps = size(visit_matrix,2);
    sample_logic = sum(visit_matrix,1);
    
    hexa_model.visits = zeros(size(visit_matrix));
    hexa_model.rewards = zeros(size(visit_matrix));
    
    p_reward = zeros(size(visit_matrix));
    p_reward(:,1) = prior(:,1);
    p_reward(:,2) = prior(:,1);
    p_stay = zeros(size(visit_matrix));
    p_stay(:,1) = prior(:,2);
    
    hexa_model.stay_go = zeros(1,size(visit_matrix,2));
    
    checked_port = 0;
    last_checked_port = checked_port;
    port_array = 1:6;

    reward_available = zeros(size(visit_matrix));
    reward_available(:,1) = 1;
    yes_reward=0;

    if sample_logic(1)==1
        checked_port = randperm(6,1);
        hexa_model.visits(checked_port,1) = 1;
        hexa_model.rewards(checked_port,1) = 1;
    end

    for t=2:max_tsteps-1

        reward_available(reward_available(:,t)==0,t) = rew_sched(reward_available(:,t)==0,t);
        reward_available(:,t+1) = reward_available(:,t);

        p_reward(:,t)   = p_reward(:,t-1);
        p_stay(:,t)     = p_stay(:,t-1);

        
       % should we check any port at this time point
       if sample_logic(t)==1
           
            % find out what ports have a reward available
            avail_array = find(reward_available(:,t)==1);

            if numel(avail_array)>=1
                % choose amongst those
                % [~,checked_port] = max(portq(avail_array));
                checked_port = randsample(avail_array,1);
            else
                % if none choose randomly            
                checked_port = randsample(port_array,1);
            end

            if last_checked_port>0  
                if checked_port == last_checked_port
                    hexa_model.stay_go(t) = 1;
                else
                    hexa_model.stay_go(t) = 0;                    
                end
            end

            hexa_model.visits(checked_port,t) = 1;   
      
           % Was the check rewarded?
           if reward_available(checked_port,t)==1
               hexa_model.rewards(checked_port,t) = 1;
               reward_available(checked_port,t+1) = 0;
               yes_reward = 1;           
           else
               yes_reward = 0;
           end
    
           last_checked_port = checked_port;

       end
       
    end

tmp                     = find(sum(visit_matrix,1)==1);
[~,visit_list_data]     = max(visit_matrix(:,tmp),[],1);    
[trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data(1:end),0,1);
[~,visit_list_model]    = max(hexa_model.visits(:,tmp),[],1);
[trans_mat_model]       = HX_ComputeTransitionMatrix(visit_list_model(1:end),0,1);
[trans_mat_rand]        = HX_ComputeTransitionMatrix( randsample(1:6, numel(visit_list_data), true, ones(1,6)), 0, 1);

exag_map = TNC_CreateRBColormap(8,'exag');
figure(249); clf; 
subplot(131); imagesc(trans_mat_data,[0 0.25]); colormap(exag_map);
subplot(132); imagesc(trans_mat_model,[0 0.25]); title('model');
subplot(133); imagesc(trans_mat_rand,[0 0.25]); title('model');

%--------
% Maybe just update the fit to -LL of observed choices?
% I think that would just be log( observed visit matrix - estimated
% probability of choices from simulations ) -> summed over total visits or
% mean per visit

trans_r2                = corr2(trans_mat_data,trans_mat_model);
trans_r2_rand           = corr2(trans_mat_data,trans_mat_rand);

all_visits              = find(sample_logic==1);
rew_logic               = sum(hexa_model.rewards,1);
all_rewards             = rew_logic(all_visits);
income_model            = movmean(all_rewards,51);
rho                     = sqrt( mean( (income_model-income).^2 ) );

income_r2               = rho;
figure(250); clf; plot(income); hold on; plot(income_model); 
title(['RMSE: ' num2str(rho)]); axis([0 numel(income_model) 0 1]); box off;

visits_for_LL = hexa_model.visits(:,all_visits);
rewards_for_LL = hexa_model.rewards(:,all_visits);

drawnow;
