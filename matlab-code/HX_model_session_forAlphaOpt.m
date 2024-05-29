function [trans_r2, income_r2] = HX_model_session_forAlphaOpt(x1,x2,x3,x4,x5)
% Creating a simplified version of model code to allow optimization of
% alpha as a function of tau1 and tau2

    global visit_matrix
    global cost_per_port
    global rew_sched
    global income

    epsilon = 0.1;

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
    p_reward(:,1:2) = 1/6;
    p_stay = zeros(size(visit_matrix));
    p_stay(:,1) = 1/6;
    
    hexa_model.stay_go = zeros(1,size(visit_matrix,2));
    
    checked_port = 0;
    last_checked_port = checked_port;
    port_array = 1:6;

    reward_available = zeros(size(visit_matrix));
    reward_available(:,1) = 1;
    yes_reward=0;

    v_ind = 1:sum(sample_logic);
    alpha_vis = x1 + (x2*(1-exp(-v_ind/x4)) .* (x3*exp(-v_ind/x5)));


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
           
          vis_cnt = sum(sample_logic(1,1:t),2);
    
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
            p_stay(checked_port,t)     = alpha_vis(vis_cnt)*yes_reward + (1-alpha_vis(vis_cnt))*p_stay(checked_port,t-1);
        else
            p_reward(checked_port,t)   = alpha_vis(vis_cnt)*yes_reward + (1-alpha_vis(vis_cnt))*p_reward(checked_port,t-1);
        end

       end
           last_checked_port = checked_port;
    end

tmp                     = find(sum(visit_matrix,1)==1);
[~,visit_list_data]     = max(visit_matrix(:,tmp),[],1);    
[trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data,0,1);
[~,visit_list_model]    = max(hexa_model.visits(:,tmp),[],1);    
[trans_mat_model]       = HX_ComputeTransitionMatrix(visit_list_model,0,1);

trans_r2                = corr2(trans_mat_data,trans_mat_model);

all_visits              = find(sample_logic==1);
rew_logic               = sum(hexa_model.rewards,1);
all_rewards             = rew_logic(all_visits);
income_model            = cumsum(all_rewards);


rho                     = sqrt( mean( (income_model-income).^2 ) );

income_r2               = rho;
figure(250); clf; plot(income); hold on; plot(income_model); 
title(['RMSE: ' num2str(rho)]); axis([0 numel(income_model) 0 numel(income_model)/2]); box off;
drawnow;
