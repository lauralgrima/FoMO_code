function [hexa_model] = HX_model_session(hexa_data_an,plot_out)

hexa_model.seed = randperm(100,1);
rng(hexa_model.seed);

% sampling rate is now set to camera frame rate
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
policy.type = 'e-greedy'; % out of type = {'softmax','greedy','e-greedy','random','proportional','e-proportional'}
policy.params.epsilon = 0.075;

belief.type = 'win-stay'; % out of type = {'win-stay','proportional','kernel','spatial','pdf','pdf-space'}
% 'win-stay' - biased towards staying at current port after reward; visit with no reward explores
% 'matching' - P(rew|port) = sum(rew(port))
% 'match-shift' - P(rew|port) = sum(rew(port)) +
%           tendency to shift after a success
% 'matchP-shift' - P(rew|port) = num_rew./num_visits +
%           tendency to shift after a success
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
reward_available = zeros(size(hexa_data_an.visits));
reward_availableI = zeros(size(hexa_data_an.visits));
hexa_model.ideal = zeros(1,max_tsteps);
p_reward = zeros(size(hexa_data_an.visits));
p_reward(:,1) = 1/6;
p_NOreward = zeros(size(hexa_data_an.visits));


for t=2:max_tsteps-1
    
   reward_available(reward_available(:,t)==0,t) = hexa_model.rew_sched(reward_available(:,t)==0,t);
   reward_available(:,t+1) = reward_available(:,t);

   p_reward(:,t) = p_reward(:,t-1);

   % should we check any port at this time point
   if sample_logic(t)==1

       if max(reward_available(:,t))==1
           hexa_model.ideal(t) = 1;
       end
       
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
                    p_reward(checked_port,t) = p_reward(checked_port,t).*0.67;
                end

           case 'matchP-shift' %- P(rew|port) = sum(rew(port))./sum(visits)               
               p_reward(:,t) = (sum(hexa_model.rewards(:,1:t),2)+0.16) ./ (sum(hexa_model.visits(:,1:t),2)+1);
                if yes_reward
                    p_reward(checked_port,t) = p_reward(checked_port,t).*0.67;
                end

           case 'matchP-shift-local' %- P(rew|port) = sum(rew(port))./sum(visits)
               p_reward(checked_port,t) = (sum(hexa_model.rewards(checked_port,1:t),2)+0.16) ./ (sum(hexa_model.visits(checked_port,1:t),2)+1);
                if yes_reward
                    p_reward(checked_port,t) = p_reward(checked_port,t).*0.67;
                end

           case 'matchP-shift-spatial' %- proportional + discount due to distance to port from current location
               p_reward(:,t) = (sum(hexa_model.rewards(:,1:t),2)+0.16) ./ (sum(hexa_model.visits(:,1:t),2)+1);
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
    
       % Compute the idealized choice model
       reward_availableI(reward_availableI(:,t)==0,t) = hexa_model.rew_sched(reward_availableI(:,t)==0,t);
       reward_availableI(:,t+1) = reward_availableI(:,t);

       % should we check any port at this time point
       if sample_logic(t)==1

           tmp = find(reward_availableI(:,t)==1);
           if numel(tmp)>0
               hexa_model.ideal(t) = 1;
               checked_port = tmp(randperm(numel(tmp),1));

               % Was the check rewarded?
               if reward_availableI(checked_port,t)==1
                   reward_availableI(:,t+1) = reward_availableI(:,t);
                   reward_availableI(checked_port,t+1) = 0;
               end       
           end

       end

   end
    
end

% -----------------------------------------------------------------
% PLOTTING OUTPUT GRAPHS AND SUMMARY STATS FROM MODEL
% -----------------------------------------------------------------
disp(['Model rewards collected: ' num2str(sum(sum(hexa_model.rewards))) ' ; ' num2str(100*(sum(sum(hexa_model.rewards)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Mouse rewards collected: ' num2str(sum(sum(hexa_data_an.rewards))) ' ; ' num2str(100*(sum(sum(hexa_data_an.rewards)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Idealized Mouse rewards: ' num2str(sum(hexa_model.ideal)) ' ; ' num2str(100*(sum(hexa_model.ideal))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Max rewards available: ' num2str(sum(sum(hexa_model.rew_sched)))])

if plot_out
    figure(60); imagesc(p_reward);
    
    hexa_model.port_chk_t = find(sum(hexa_model.visits,1)==1);
    hexa_model.port_chk   = hexa_model.visits(:,hexa_model.port_chk_t);
    
    hexa_model.trans_mat = zeros(6,6);
    
    % compute the transition matrix of visits
    for qq=1:6
        tmp = find(hexa_model.port_chk(qq,1:end-1)==1);
        trans_cnt = zeros(1,6);
        for pp=tmp
            next_visit = find(hexa_model.port_chk(:,pp+1)==1);
            trans_cnt(next_visit) = trans_cnt(next_visit) + 1;
        end
        hexa_model.trans_mat(qq,:) = trans_cnt;
    %     hexa_model.trans_mat(:,qq) = trans_cnt';
    end
    
    % Want to compare these two matrices + reward fraction with models 
    figure(32); clf;
    subplot(1,4,1:2);
    imagesc(hexa_model.port_chk);
    xlabel('Unique visits'); ylabel('Port');
    subplot(1,4,3);
    plot(sum(hexa_model.port_chk,2),1:6,'o-'); hold on;
    plot(sum(hexa_data_an.port_chk,2),1:6,'ko-');
    legend({'Model','Data'});
    xlabel('Total visits'); ylabel('Port');
    set(gca,'YDir','reverse');
    axis([0 size(hexa_model.port_chk,2)/2 1 6]); box off;
    subplot(1,4,4);
    imagesc(hexa_model.trans_mat./size(hexa_model.port_chk,2),[0 0.15]); colormap(tmap); colorbar;
    xlabel('Port'); ylabel('Port');
    title('Transition probability');
end