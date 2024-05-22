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
p_option(1,:) = [0.5*ones(1,200) 0.75*ones(1,500) 0.25*ones(1,190) 0.66*ones(1,500) 0.33*ones(1,350)];
p_option(2,:) = [0.5*ones(1,200) 0.25*ones(1,500) 0.75*ones(1,190) 0.33*ones(1,500) 0.66*ones(1,350)];

[qual_map] = TNC_CreateRBColormap(8,'cat2');
qual_map = qual_map([1 3 5],:);

num_trials = size(p_option,2);

figure(400); clf;
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

    if jj<num_trials
        p_rew(1,jj+1) = alpha*rewards(1,jj) + beta*p_rew(1,jj);
        p_rew(2,jj+1) = alpha*rewards(2,jj) + beta*p_rew(2,jj);
    end

    
end

plot(1:num_trials,conv(visits(1,:),[0 ones(1,25) 0]/25,'same'),'color',qual_map(1,:),'LineWidth',2);
plot(1:num_trials,conv(visits(2,:),[0 ones(1,25) 0]/25,'same'),'color',qual_map(2,:),'LineWidth',2);
