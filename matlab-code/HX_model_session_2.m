function [hexa_model] = HX_model_session_2(hexa_data_an,policy_type,belief_type,cost_per_port,port_intervals,dynamic_epsilon,base_stay,plot_out)

hexa_model.seed = randperm(1000,1);
rng(hexa_model.seed);

% Interport distance matrix
hexa_model.interportdist =      ...
[0	14	18	70	72.2	65.5 ;  ...
14	0	22.8	56	65.5	42; ...
18	22.8	0	72.2	70	56; ...
70	56	72.2	0	18	22.8;   ...
72.2	65.5	70	18	0	14; ...
65.5	42	56	22.8	14	0];

% sampling rate is now set to 1 Hz
frame_rate = 1;

% Make the BigBoi reward environment
hexa_model.rew_sched = zeros(size(hexa_data_an.visits));
for qq=1:6
    hexa_model.rew_sched(qq,1:port_intervals(qq)*frame_rate:end) = 1;
end
hexa_model.rew_sched(:,2) = 1;

max_reward = sum(sum(hexa_model.rew_sched));

% Pass in data file and policy choice
policy.type = policy_type; % out of type = {'softmax','greedy','e-greedy','random','proportional','e-proportional'}
if dynamic_epsilon==1
    policy.params.epsilon = 1;
    epsilon_tau = 33;
    stable_epsilon = 0.1;
    hexa_model.epsilon_tau = epsilon_tau;
    hexa_model.stable_epsilon = stable_epsilon;
elseif dynamic_epsilon==-1
    policy.params.epsilon = 0;
    stable_epsilon = 0;
    hexa_model.epsilon_tau = 0;
    hexa_model.stable_epsilon = stable_epsilon;
else
    policy.params.epsilon = 0.1;
    stable_epsilon = 0.1;
    hexa_model.epsilon_tau = 0;
    hexa_model.stable_epsilon = stable_epsilon;
end

if strmatch(policy.type,'softmax')
    stable_epsilon = 0;
    hexa_model.stable_epsilon = stable_epsilon;
end

belief.type = belief_type; % out of type = {'win-stay','proportional','kernel','spatial','pdf','pdf-space'}
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
% 'hazard-v' - attempt to estimate true hazard in units of visits P(rew|port,visit)
% 'pdf-space' - combined belief about posterior and discounting by distance

belief.params = 0;

% Format the simulation to produce same size output
% Plot and compare model to data

max_tsteps = size(hexa_data_an.visits,2);
sample_logic = sum(hexa_data_an.visits,1);

hexa_model.visits = zeros(size(hexa_data_an.visits));
hexa_model.rewards = zeros(size(hexa_data_an.visits));
hexa_model.ideal = zeros(size(hexa_data_an.visits));
hexa_model.random = zeros(size(hexa_data_an.visits));

reward_available = zeros(size(hexa_data_an.visits));
reward_available(:,1) = 1;
reward_availableI = zeros(size(hexa_data_an.visits));
reward_availableR = zeros(size(hexa_data_an.visits));

p_reward = zeros(size(hexa_data_an.visits));
p_reward(:,1:2) = 1/6;
p_NOreward = zeros(size(hexa_data_an.visits));
p_stay = zeros(size(hexa_data_an.visits));
p_stay(:,1) = 1/6;

hexa_model.stay_go = zeros(1,size(hexa_data_an.visits,2));

checked_port = 0;
last_checked_port = checked_port;
port_array = 1:6;

for t=2:max_tsteps-1
    
    reward_available(reward_available(:,t)==0,t) = hexa_model.rew_sched(reward_available(:,t)==0,t);
    reward_available(:,t+1) = reward_available(:,t);

    reward_availableI(reward_availableI(:,t)==0,t) = hexa_model.rew_sched(reward_availableI(:,t)==0,t);
    reward_availableI(:,t+1) = reward_availableI(:,t);

    reward_availableR(reward_availableR(:,t)==0,t) = hexa_model.rew_sched(reward_availableR(:,t)==0,t);
    reward_availableR(:,t+1) = reward_availableR(:,t);
    
    p_reward(:,t)   = p_reward(:,t-1);
    p_stay(:,t)     = p_stay(:,t-1);

   % should we check any port at this time point
   if sample_logic(t)==1
       
      % stay or shift decision

      if last_checked_port>0  
          
          if rand(1) < p_stay(last_checked_port,t)

              checked_port = last_checked_port;
              % disp('Stayed...'); 
              hexa_model.stay_go(t) = 1;
              hexa_model.visits(checked_port,t) = 1;

          else % switch and use matching 
    
              hexa_model.stay_go(t) = 0;

              % Use 'policy' to govern port choice
               switch policy.type
               
                   case 'softmax'
                        if rand(1)>policy.params.epsilon
                            if max(p_reward(port_array~=last_checked_port,t))<=0
                                port_array(port_array~=last_checked_port)
                                cost_per_port(port_array~=last_checked_port,checked_port)
                                ones(1,5)./cost_per_port(port_array~=last_checked_port,checked_port)
                                checked_port = randsample(port_array(port_array~=last_checked_port),1,true,ones(5,1)./cost_per_port(port_array~=last_checked_port,checked_port));
                            else
                                tmp_p_reward = p_reward(port_array~=last_checked_port,t);
                                this_p_reward = exp(tmp_p_reward)./sum(exp(tmp_p_reward));
                                checked_port = randsample(port_array(port_array~=last_checked_port),1,true,this_p_reward./cost_per_port(port_array~=last_checked_port,checked_port));
                            end
                        else
                            checked_port = randsample(port_array(port_array~=last_checked_port),1);
                        end
                        hexa_model.visits(checked_port,t) = 1;
        
                   case 'e-proportional'   
                        if rand(1)>policy.params.epsilon
                            if max(p_reward(port_array~=last_checked_port,t))<=0
                                checked_port = randsample(port_array(port_array~=last_checked_port),1,true,[0.2 0.2 0.2 0.2 0.2]./cost_per_port(port_array~=last_checked_port,checked_port));
                            else
                                checked_port = randsample(port_array(port_array~=last_checked_port),1,true,p_reward(port_array~=last_checked_port,t)./cost_per_port(port_array~=last_checked_port,checked_port));
                            end
                        else
                            checked_port = randsample(port_array(port_array~=last_checked_port),1);
                        end
                        hexa_model.visits(checked_port,t) = 1;

                   case 'e-greedy'   
                        if rand(1)>policy.params.epsilon
                            if max(p_reward(port_array~=last_checked_port,t))<=0
                                checked_port = randsample(port_array(port_array~=last_checked_port),1,true,[0.2 0.2 0.2 0.2 0.2]./cost_per_port(port_array~=last_checked_port,checked_port));
                            else
                                [~,max_port] = max( p_reward(port_array~=last_checked_port,t)./cost_per_port(port_array~=last_checked_port,checked_port) );
                                valid_ports = port_array(port_array~=last_checked_port);
                                checked_port = valid_ports(max_port);
                            end
                        else
                            checked_port = randsample(port_array(port_array~=last_checked_port),1);
                        end
                        hexa_model.visits(checked_port,t) = 1;
                        
                   case 'random' % random policy
                       checked_port = randperm(6,1);
                       hexa_model.visits(checked_port,t) = 1;
    
               end
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
           
            if dynamic_epsilon==1
                policy.params.epsilon = stable_epsilon+exp(-sum(sum(hexa_model.rewards(:,1:t)))./epsilon_tau);
            end
           
       else
           yes_reward = 0;
       end       
       
      % Compute the omniscient choice model
       tmp = find(reward_availableI(:,t)==1);
       if numel(tmp)>0
           checked_portI = tmp(randperm(numel(tmp),1));
           hexa_model.ideal(checked_portI,t) = 1;
           reward_availableI(checked_portI,t+1) = 0;
       end

      % Compute the random choice model
      checked_portR = randperm(6,1);
       if reward_availableR(checked_portR,t)==1
           hexa_model.random(checked_portR,t) = 1;
           reward_availableR(:,t+1) = reward_availableR(:,t);
           reward_availableR(checked_portR,t+1) = 0;
       end       


       % Update belief { Pr(R|port,t) } according to different models
       switch belief.type

           case 'p_check_match' %- attempt to estimate P(rew|port,t)

               % p_check part
               for port_id=1:6
                   if numel(find(hexa_model.visits(port_id,1:t)==1 & hexa_model.stay_go(1:t)==1))==0
                       stay_rewards = 0;
                   else
                    stay_rewards = sum(hexa_model.rewards(port_id,hexa_model.visits(port_id,1:t)==1 & hexa_model.stay_go(1:t)==1)) ./ numel(find(hexa_model.visits(port_id,1:t)==1 & hexa_model.stay_go(1:t)==1));
                   end
                   if numel(find(hexa_model.visits(port_id,1:t)==1 & hexa_model.stay_go(1:t)==0)) == 0
                       return_rewards = 0;
                   else
                    return_rewards = sum(hexa_model.rewards(port_id,hexa_model.visits(port_id,1:t)==1 & hexa_model.stay_go(1:t)==0)) ./ numel(find(hexa_model.visits(port_id,1:t)==1 & hexa_model.stay_go(1:t)==0));
                   end

                   p_reward(port_id,t)   = return_rewards;  
                   p_stay(port_id,t)     = stay_rewards;

               end

                % should switch to 
                

                if yes_reward
                    p_stay(checked_port,t) = base_stay;
                                % NOTE for N-2 transition to work would need to also modulate
                                % p_reward, not just p_stay
                else
                    % nothing...
                    if p_stay(checked_port,t) < base_stay
                        p_stay(checked_port,t) = base_stay;
                    end
                end



                % % p_stay part
                % for qq=1:6
                %     if numel(find(hexa_model.rewards(qq,1:t)==1,1,'last'))>0
                %         all_rew = find(hexa_model.rewards(qq,1:t)==1);
                %         if numel(all_rew)>1
                %             hazard = cumsum(diff(all_rew))./sum(diff(all_rew)); % choose tau such that by ~5 tau full expectation is returned
                %         else
                %             hazard = 2;
                %         end
                %         since_last       = t-find(hexa_model.rewards(qq,1:t)==1,1,'last')
                % 
                %         % if since_last > numel(hazard) | since_last==0
                %             p_stay(qq,t)  = p_reward(qq,t);
                %         % else
                %         %     p_stay(qq,t)  = p_reward(qq,t) .* hazard(since_last);
                %         % end
                %     else
                %         p_stay(qq,t)  = p_reward(qq,t);
                %     end
                % end

            case 'WSLS' %- attempt to estimate P(rew|port,t)
                if yes_reward
                    p_stay(:,t) = policy.params.epsilon;
                    p_stay(checked_port,t) = 1-policy.params.epsilon;
                    p_reward(:,t) = 0.16;
                else
                    % nothing...
                    p_stay(:,t) = policy.params.epsilon;
                    p_reward(:,t) = 0.16;
                end

           otherwise % do nothing

      end    
   end
   last_checked_port = checked_port;
end

% -----------------------------------------------------------------
% PLOTTING OUTPUT GRAPHS AND SUMMARY STATS FROM MODEL
% -----------------------------------------------------------------
disp(['Model rewards collected: ' num2str(sum(sum(hexa_model.rewards))) ' ; ' num2str(100*(sum(sum(hexa_model.rewards)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Mouse rewards collected: ' num2str(sum(sum(hexa_data_an.rewards))) ' ; ' num2str(100*(sum(sum(hexa_data_an.rewards)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Omniscient policy rewards: ' num2str(sum(sum(hexa_model.ideal))) ' ; ' num2str(100*(sum(sum(hexa_model.ideal)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Random policy rewards: ' num2str(sum(sum(hexa_model.random))) ' ; ' num2str(100*(sum(sum(hexa_model.random)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Max rewards available: ' num2str(sum(sum(hexa_model.rew_sched)))])

if plot_out
    
    pfm = figure(61); clf; subplot(131);
    % plot income over time slope vs omniscient and random
    plot(cumsum(sum(hexa_model.ideal,1)),'k-','linewidth',2); hold on;
    plot(cumsum(sum(hexa_model.random,1)),'k-','linewidth',1); hold on;
    plot(cumsum(sum(hexa_data_an.rewards,1)),'-','color',[1 0 0.33],'linewidth',2); hold on;
    plot(cumsum(sum(hexa_model.rewards,1)),'-','color',[0.5 0 0.16],'linewidth',2); hold on;
%     legend({'Omniscient','Random',['Data: ' hexa_data_an.filename(1:end-3) ' ; s' num2str(hexa_data_an.session)] ,['Model: ' belief.type ' & ' policy.type]},'Location','Northwest');
    legend({'Omniscient','Random','Data','Model'},'Location','Northwest');
    box off;
    ylabel('Cumulative rewards'); xlabel('Session Time');
    axis([0 size(hexa_model.ideal,2) 0 max(cumsum(sum(hexa_model.ideal,1)))]);

    
    % compute sliding window slope estimate
    ideal_income = cumsum(sum(hexa_model.ideal,1));
    mouse_income = cumsum(sum(hexa_data_an.rewards,1));
    model_income = cumsum(sum(hexa_model.rewards,1));
    random_income = cumsum(sum(hexa_model.random,1));
    clear hexa_model.slope;
    cnt=1;
    for ww=2000:50:numel(ideal_income)
        xs = ww-1999:ww;
        est = polyfit(xs,ideal_income(xs),1);
        hexa_model.slope.ideal(cnt) = est(1);
        
        est = polyfit(xs,mouse_income(xs),1);
        hexa_model.slope.mouse(cnt) = est(1);
        
        est = polyfit(xs,model_income(xs),1);
        hexa_model.slope.model(cnt) = est(1);
        
        est = polyfit(xs,random_income(xs),1);
        hexa_model.slope.random(cnt) = est(1);
        
        hexa_model.slope.x(cnt) = ww;
        
        cnt = cnt+1;
    end
    subplot(132);
    plot(hexa_model.slope.x,hexa_model.slope.ideal,'k-','linewidth',2); hold on;
    plot(hexa_model.slope.x,hexa_model.slope.random,'k-','linewidth',1); hold on;
    plot(hexa_model.slope.x,hexa_model.slope.mouse,'-','color',[1 0 0.33],'linewidth',2); hold on;
    plot(hexa_model.slope.x,hexa_model.slope.model,'-','color',[0.5 0 0.16],'linewidth',2); hold on;
%     legend({'Omniscient','Random',['Data: ' hexa_data_an.filename(1:end-3) ' ; s' num2str(hexa_data_an.session)] ,['Model: ' belief.type ' & ' policy.type]},'Location','Northwest');
    box off;
    ylabel('Local Income'); xlabel('Session Time');
%     axis([0 size(hexa_model.ideal,2) 0.5*max(cumsum(sum(hexa_model.random,1)))./size(hexa_model.random,2) 1.25*max(cumsum(sum(hexa_model.ideal,1)))./size(hexa_model.ideal,2)]);
    axis([0 size(hexa_model.ideal,2) 0.015 0.075]);

    subplot(133);
    plot(1:6,sum(hexa_data_an.visits,2),'o-','color',[1 0 0.33],'linewidth',2); hold on;
    plot(1:6,sum(hexa_model.visits,2),'o-','color',[0.5 0 0.16],'linewidth',2);
%     legend({['Data: ' hexa_data_an.filename(1:end-3) ' ; s' num2str(hexa_data_an.session)] ,['Model: ' belief.type ' & ' policy.type]},'Location','Northwest');
    xlabel('Total visits'); ylabel('Port');
    axis([0 7 0 max(sum(hexa_data_an.visits,2))*1.5]); 
    box off;

    orient(pfm,'landscape');
    print(pfm, ['~/Data_' hexa_data_an.filename(1:end-3) '_s' num2str(hexa_data_an.session) '_Model_' belief.type '_' policy.type], '-dpdf','-bestfit');

end

hexa_model.p_reward = p_reward;
hexa_model.p_stay = p_stay;
hexa_model.stable_epsilon = stable_epsilon;
