function [hexa_model] = HX_model_session(hexa_data_an,policy_type,belief_type,plot_out)

hexa_model.seed = randperm(100,1);
rng(hexa_model.seed);

% Interport distance matrix
hexa_model.interportdist = ...
[0	14	18	70	72.2	65.5 ;...
14	0	22.8	56	65.5	42;...
18	22.8	0	72.2	70	56;...
70	56	72.2	0	18	22.8;...
72.2	65.5	70	18	0	14;...
65.5	42	56	22.8	14	0];

% sampling rate is now set to 1 Hz
frame_rate = 1;

% Make the BigBoi reward environment
hexa_model.rew_sched = zeros(size(hexa_data_an.visits));
hexa_model.rew_sched(1,1:30*frame_rate:end) = 1;
hexa_model.rew_sched(2,200:60*frame_rate:end) = 1;
hexa_model.rew_sched(3,400:240*frame_rate:end) = 1;
hexa_model.rew_sched(4,600:240*frame_rate:end) = 1;
hexa_model.rew_sched(5,800:1200*frame_rate:end) = 1;
hexa_model.rew_sched(6,1000:2400*frame_rate:end) = 1;
hexa_model.rew_sched(:,1) = 1;

max_reward = sum(sum(hexa_model.rew_sched));

% Pass in data file and policy choice
policy.type = policy_type; % out of type = {'softmax','greedy','e-greedy','random','proportional','e-proportional'}
policy.params.epsilon = 0.075;

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
reward_availableI = zeros(size(hexa_data_an.visits));
reward_availableR = zeros(size(hexa_data_an.visits));

p_reward = zeros(size(hexa_data_an.visits));
p_reward(:,1) = 1/6;
p_NOreward = zeros(size(hexa_data_an.visits));


for t=2:max_tsteps-1
    
    reward_available(reward_available(:,t)==0,t) = hexa_model.rew_sched(reward_available(:,t)==0,t);
    reward_available(:,t+1) = reward_available(:,t);

    reward_availableI(reward_availableI(:,t)==0,t) = hexa_model.rew_sched(reward_availableI(:,t)==0,t);
    reward_availableI(:,t+1) = reward_availableI(:,t);

    reward_availableR(reward_availableR(:,t)==0,t) = hexa_model.rew_sched(reward_availableR(:,t)==0,t);
    reward_availableR(:,t+1) = reward_availableR(:,t);
    
    p_reward(:,t) = p_reward(:,t-1);

   % should we check any port at this time point
   if sample_logic(t)==1
       
      % Use 'policy' to govern port choice
       switch policy.type
           
           case 'softmax'
                p_choice = softmax(p_reward(:,t)./sum(p_reward(:,t)));
                checked_port = find(rand(1)>cumsum(p_choice),1,'last')+1;
                if numel(checked_port)==0
                    checked_port=1;
                end
                hexa_model.visits(checked_port,t) = 1;

           case 'greedy'   
                [~,checked_port] = max(p_reward(:,t));
                hexa_model.visits(checked_port,t) = 1;
          
           case 'e-greedy'   
                if rand(1)>policy.params.epsilon
                    [~,checked_port] = max(p_reward(:,t));
                else
                    checked_port = randperm(6,1);
                end
                hexa_model.visits(checked_port,t) = 1;

           case 'proportional'   
                checked_port = randsample(1:6,1,true,p_reward(:,t));
                hexa_model.visits(checked_port,t) = 1;

           case 'e-proportional'   
                if rand(1)>policy.params.epsilon
                    checked_port = randsample(1:6,1,true,p_reward(:,t));
                else
                    checked_port = randperm(6,1);
                end
                hexa_model.visits(checked_port,t) = 1;
                
           case 'random' % random policy
               checked_port = randperm(6,1);
               hexa_model.visits(checked_port,t) = 1;

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
           case 'win-stay' %- biased towards staying at current port after reward; visit with no reward explores
                p_reward(:,t) = 0.02;
                if yes_reward
                    p_reward(checked_port,t) = 0.9;
                else
                    p_reward(checked_port,t) = 0.02;
                end

           case 'matching' %- P(rew|port) = sum(rew(port))
               p_reward(:,t) = sum(hexa_model.rewards(:,1:t),2)+1;

           case 'match-shift' %- P(rew|port) = sum(rew(port))
               p_reward(:,t) = sum(hexa_model.rewards(:,1:t),2)+1;
                if yes_reward
                    p_reward(checked_port,t) = 1/300;
                end

           case 'matchP-shift' %- P(rew|port) = sum(rew(port))./sum(visits)               
               p_reward(:,t) = (sum(hexa_model.rewards(:,1:t),2)+0.16) ./ (sum(hexa_model.visits(:,1:t),2)+1);
                if yes_reward
                    p_reward(checked_port,t) = 1/300;
                end

           case 'matchP-shift-local' %- P(rew|port) = sum(rew(port))./sum(visits)
               p_reward(checked_port,t) = (sum(hexa_model.rewards(checked_port,1:t),2)+0.16) ./ (sum(hexa_model.visits(checked_port,1:t),2)+1);
                if yes_reward
                    p_reward(checked_port,t) = 1/300;
                end

           case 'matchP-shift-spatial' %- proportional + discount due to distance to port from current location
               p_reward(:,t) = (sum(hexa_model.rewards(:,1:t),2)+0.16) ./ (sum(hexa_model.visits(:,1:t),2)+1);
                if yes_reward
                    p_reward(:,t) = p_reward(:,t) ./ hexa_model.interportdist(:,checked_port);
                    p_reward(checked_port,t) = 1/300;
                end
                
           case 'matchP-shift-local-spatial' %- P(rew|port) = sum(rew(port))./sum(visits)
               p_reward(checked_port,t) = (sum(hexa_model.rewards(checked_port,1:t),2)+0.16) ./ (sum(hexa_model.visits(checked_port,1:t),2)+1);
                if yes_reward
                    p_reward(:,t) = p_reward(:,t) ./ hexa_model.interportdist(:,checked_port);
                    p_reward(checked_port,t) = 1/300;
                end
                

           case 'hazard' %- attempt to estimate P(rew|port,t)

                for zzz=1:6
                                   
                        all_rew_t = find(hexa_data_an.rewards(zzz,1:t)==1);
                        all_vis_t = find(hexa_data_an.visits(zzz,1:t)==1);
                
                        if numel(all_rew_t)>0
                        
                            if numel(find(t==all_rew_t))==1
                                pdf                = hist([all_rew_t(1) diff(all_rew_t)],0:1:1e6);
                                cdf                = cumsum(pdf)./sum(pdf);                
                            end
                
                            if t-all_rew_t(end) < 1e6
                                p_reward(zzz,t) = cdf(t-all_rew_t(end)+1); % i.e. the cdf evaluated at current t
                            else
                                p_reward(zzz,t) = cdf(end);
                            end
                
                        elseif numel(find(t==all_vis_t))>0 % also get the negative cdf (i.e. visit intervals that fail to elicit reward)
                            
                            pdf                = hist([all_vis_t(1) diff(all_vis_t)],0:1:1e6);
                            cdf                = cumsum(pdf)./sum(pdf);
                
                            if t-all_vis_t(end) < 1e6
                                p_NOreward(zzz,t) = cdf(t-all_vis_t(end)+1); % i.e. the cdf evaluated at current t
                            else
                                p_NOreward(zzz,t) = cdf(end);
                            end
                        
                        end
                
                end


           otherwise % do nothing

       end
    
   end
    
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
    legend({'Omniscient','Random',['Data: ' hexa_data_an.filename(1:end-3) ' ; s' num2str(hexa_data_an.session)] ,['Model: ' belief.type ' & ' policy.type]},'Location','Northwest');
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
    legend({'Omniscient','Random',['Data: ' hexa_data_an.filename(1:end-3) ' ; s' num2str(hexa_data_an.session)] ,['Model: ' belief.type ' & ' policy.type]},'Location','Northwest');
    box off;
    ylabel('Local Income'); xlabel('Session Time');
    axis([0 size(hexa_model.ideal,2) 0.5*max(cumsum(sum(hexa_model.ideal,1)))./size(hexa_model.ideal,2) 1.25*max(cumsum(sum(hexa_model.ideal,1)))./size(hexa_model.ideal,2)]);

    subplot(133);
    plot(1:6,sum(hexa_data_an.visits(:,round(end/2):end),2),'o-','color',[1 0 0.33],'linewidth',2); hold on;
    plot(1:6,sum(hexa_model.visits(:,round(end/2):end),2),'o-','color',[0.5 0 0.16],'linewidth',2);
    legend({['Data: ' hexa_data_an.filename(1:end-3) ' ; s' num2str(hexa_data_an.session)] ,['Model: ' belief.type ' & ' policy.type]},'Location','Northwest');
    xlabel('Total visits (2nd half)'); ylabel('Port');
    axis([0 7 0 max(sum(hexa_data_an.visits(:,end/2:end),2))*1.5]); 
    box off;

    print(pfm, ['Data_' hexa_data_an.filename(1:end-3) '_s' num2str(hexa_data_an.session) '_Model_' belief.type '_' policy.type], '-dpdf', '-bestfit', '-vector');

end