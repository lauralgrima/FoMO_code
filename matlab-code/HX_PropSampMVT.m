%% Script to examine what MVT schedules look like under the ALLOCATE model 

% Laura version (not used here at the moment)
% Environment: initial mean forage time * (depletion rate ^ (reward number in patch-1)).

% Construct a depleting patch environment 

[qual_map] = TNC_CreateRBColormap(8,'cat2');
qual_map = qual_map([1 3 5],:);
figure(300); clf;

P_rew_tau = 8;

for jj=1:3
    P_rew_init_poss = [1 0.75 0.5];
    P_rew_init = P_rew_init_poss(jj);
    rewards = 0:20;    
    P_rew = P_rew_init * exp(-rewards/P_rew_tau);
    
    subplot(141);
    plot(rewards,P_rew,'color',qual_map(jj,:)); hold on; box off;
    scatter(rewards,P_rew,50,qual_map(jj,:),'filled');
    axis([-1 21 0 1]); ylabel('P(rew)'); xlabel('Visits');

end

num_agents  = 20;
num_iters   = 75;

Rew_To_Leave = zeros(num_iters,3,num_agents);
Vis_To_Leave = zeros(num_iters,3,num_agents);

for qq=1:num_agents

P_leave = 0.05;

P_leave     = 0.05 + randn(1)./50;
if P_leave<0
    P_leave = 1e-3;
end

    for jj=1:3
    
        P_rew_init = P_rew_init_poss(jj);
    
        for iter=1:num_iters
            
            rew_cnt = 0;
            visit   = 1;
            P_stay  = [];
            stay    = 1;
    
            while stay==1
    
                P_rew = P_rew_init * exp(-rew_cnt/P_rew_tau);
    
                % visit rewarded?
                if rand(1)<P_rew
                    rew_cnt=rew_cnt+1;
                end
    
                % stay or leave?
                P_stay(visit) = (1+rew_cnt) ./ visit;
    
                stay = randsample([1 2],1,true,[P_stay(visit) P_leave]);
                % stay = rand(1) < P_stay(visit);
                
                if stay==1
                    visit = visit+1;
                else
                    Rew_To_Leave(iter,jj,qq) = rew_cnt;
                    Vis_To_Leave(iter,jj,qq) = visit;
                end
    
            end
            % 
            % subplot(144);
            % plot(1:visit,P_stay,'Color',qual_map(jj,:)); hold on;
            % ylabel('P(stay)'); xlabel('Visits');
            % drawnow;
        
        end
    
    end

end

subplot(142);
boxplot(squeeze(mean(Rew_To_Leave,3)));
axis([0 4 0 12]);
ylabel('Rewards before leave'); xlabel('Patch richness');
box off;

subplot(143);
boxplot(squeeze(mean(Vis_To_Leave,3)));
axis([0 4 0 28]);
ylabel('Visits before leave'); xlabel('Patch richness');
box off;

%% Examine a matching income task like Sugrue/Bari

clear p_option; 

alpha = 0.1
beta = 1-alpha;
epsilon = 0.1

% Create the environment 
p_option(1,:) = [0.5*ones(1,200) 0.9*ones(1,500) 0.1*ones(1,190) 0.7*ones(1,500) 0.3*ones(1,350) 0.55*ones(1,500) 0.45*ones(1,350)];
p_option(2,:) = [0.5*ones(1,200) 0.1*ones(1,500) 0.9*ones(1,190) 0.3*ones(1,500) 0.7*ones(1,350) 0.45*ones(1,500) 0.55*ones(1,350)];

[qual_map] = TNC_CreateRBColormap(8,'cat2');
qual_map = qual_map([1 3 5],:);

num_trials = size(p_option,2);

figure(400); clf;
subplot(1,3,1:2);
plot(1:num_trials,p_option(1,:),'color',qual_map(1,:),'LineWidth',4); hold on;
plot(1:num_trials,p_option(2,:),'color',qual_map(2,:),'LineWidth',4); box off;
axis([0 num_trials 0 1]);

p_rew = zeros(2,num_trials);
p_rew(:,1) = 0.5;

visits = zeros(2,num_trials);
rewards = zeros(2,num_trials);

for jj=1:num_trials

    if rand(1)<epsilon
        choice = randperm(2,1);
    else
        choice = randsample([1 2],1,true,[p_rew(1,jj) p_rew(2,jj)]);
    end

    visits(choice,jj) = 1;

    if rand(1)<p_option(choice,jj)
        rewards(choice,jj)=1;
    end

    p_rew(:,jj+1) = p_rew(:,jj);

    if jj<num_trials
        p_rew(choice,jj+1) = alpha*rewards(choice,jj) + beta*p_rew(choice,jj);
    end

    
end

win_ln = 33;

subplot(1,3,1:2);
plot(1:num_trials,conv(visits(1,:),[0 ones(1,win_ln) 0]/win_ln,'same'),'color',qual_map(1,:),'LineWidth',2); hold on;
plot(1:num_trials,conv(visits(2,:),[0 ones(1,win_ln) 0]/win_ln,'same'),'color',qual_map(2,:),'LineWidth',2);

subplot(1,3,3);
all_probs = unique(p_option)';
for qq=all_probs
    tmp1 = find(p_option(1,:)==qq);
    tmp2 = find(p_option(2,:)==qq);
    all_vis(1,find(all_probs==qq)) = mean(visits(1,tmp1));
    all_vis(2,find(all_probs==qq)) = mean(visits(2,tmp2));
end
plot([0 1],[0 1],'k-'); hold on;
scatter([all_probs all_probs],[all_vis(1,:) all_vis(2,:)],100,[ones(1,numel(all_probs)) 2*ones(1,numel(all_probs))],'filled'); colormap(qual_map([1 2],:));
xlabel('Option P(rew)'); ylabel('Option P(visit)');


