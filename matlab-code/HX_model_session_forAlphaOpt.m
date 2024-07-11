function [trans_r2, income_r2, visits_for_LL, rewards_for_LL] = HX_model_session_forAlphaOpt(x1,x2,x3,x4,x5,alpha_version,visit_matrix,cost_per_port,rew_sched,income,prior)
% Creating a simplified version of model code to allow optimization of
% alpha as a function of tau1 and tau2

    epsilon = 0.05;

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

    v_ind = 1:sum(sample_logic);
    switch alpha_version
        case 'doub_exp'
            alpha_vis = x1 + (x2*(1-exp(-v_ind/x4)) .* (x3*exp(-v_ind/x5)));

        case 'sig_exp'
            alpha_vis = x1 + (x2 ./ (1+exp((x4-v_ind)/(x4./6)))) .*  (x3*exp(-v_ind/x5));
            figure(9); clf; plot(alpha_vis);
    end

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
           
          rew_cnt = 1+sum( sum(hexa_model.rewards,1) , 2 );
          % rew_cnt = sum( sum(hexa_model.visits,1) , 2 );
    
          if last_checked_port>0  
              
              % stay or shift decision
              if rand(1) < p_stay(last_checked_port,t)
    
                  checked_port = last_checked_port;
                  % disp('Stayed...'); 
                  hexa_model.stay_go(t) = 1;
                  hexa_model.visits(checked_port,t) = 1;
    
              else % switch and use matching 
    
                  hexa_model.stay_go(t) = 0;
    
                  % Use 'policy' to govern port choice
                    if rand(1)>epsilon
                        if max(p_reward(port_array~=last_checked_port,t))<=0
                            checked_port = randsample(port_array(port_array~=last_checked_port),1,true,ones(5,1)./cost_per_port(port_array~=last_checked_port,checked_port));
                        else
                            checked_port = randsample(port_array(port_array~=last_checked_port),1,true,p_reward(port_array~=last_checked_port,t)./cost_per_port(port_array~=last_checked_port,checked_port));
                        end
                    else
                        checked_port = randsample(port_array(port_array~=last_checked_port),1);
                    end
                    hexa_model.visits(checked_port,t) = 1;    

              end
      
          else % first check

              checked_port = randperm(6,1);

          end

       % Was the check rewarded?
       if reward_available(checked_port,t)==1
           hexa_model.rewards(checked_port,t) = 1;
           reward_available(:,t+1) = reward_available(:,t);
           reward_available(checked_port,t+1) = 0;
           yes_reward = 1;           
       else
           yes_reward = 0;
       end

       % Update belief { Pr(R|port,t) } according to different models
       p_reward(:,t)   = p_reward(:,t-1);
       p_stay(:,t)     = p_stay(:,t-1);
        if hexa_model.stay_go(t)==1
            p_stay(checked_port,t)     = alpha_vis(rew_cnt)*yes_reward + (1-alpha_vis(rew_cnt))*p_stay(checked_port,t-1);
        else
            p_reward(checked_port,t)   = alpha_vis(rew_cnt)*yes_reward + (1-alpha_vis(rew_cnt))*p_reward(checked_port,t-1);
        end

       end
           last_checked_port = checked_port;
    end

tmp                     = find(sum(visit_matrix,1)==1);
[~,visit_list_data]     = max(visit_matrix(:,tmp),[],1);    
[trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data(1:end),0,1);
[~,visit_list_model]    = max(hexa_model.visits(:,tmp),[],1);
[trans_mat_model]       = HX_ComputeTransitionMatrix(visit_list_model(1:end),0,1);

exag_map = TNC_CreateRBColormap(8,'exag');
figure(249); clf; 
subplot(121); imagesc(trans_mat_data,[0 0.25]); colormap(exag_map);
subplot(122); imagesc(trans_mat_model,[0 0.25]); title('model');

%--------
% Maybe just update the fit to -LL of observed choices?
% I think that would just be log( observed visit matrix - estimated
% probability of choices from simulations ) -> summed over total visits or
% mean per visit

trans_r2                = corr2(trans_mat_data,trans_mat_model);

all_visits              = find(sample_logic==1);
rew_logic               = sum(hexa_model.rewards,1);
all_rewards             = rew_logic(all_visits);
income_model            = movmean(all_rewards,51);
rho                     = sqrt( mean( (income_model-income).^2 ) );
income_r2               = rho;

figure(250); clf; plot(income,'k'); hold on; plot(income_model,'r'); 
title(['RMSE: ' num2str(rho)]); axis([0 numel(income_model) 0 1]); box off;

visits_for_LL = hexa_model.visits(:,all_visits);
rewards_for_LL = hexa_model.rewards(:,all_visits);

drawnow;
