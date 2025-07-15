visit_index = 1:1:1600;
a_tau = 83;
a_scale = 0.03;
alpha_by_visit = (a_scale * exp(-(visit_index-1)./a_tau)) + 0.01;

reward_sched = [0.59,0.5,0.37,0.36,0.17,0.15];
% reward_sched = reward_sched ./ sum(reward_sched)

[port_colmap] = TNC_CreateRBColormap(6,'grima');

p_reward    = zeros(6,numel(visit_index));
p_stay      = zeros(6,numel(visit_index));
yes_reward  = zeros(1,numel(visit_index));

p_reward(:,1)   = 1/6;
p_stay(:,1)     = 0;

epsilon = 0.1;

opto_on = 0;
live_plot = 0;
decay_model = 'minus_alpha'

figure(1); clf;
% subplot(1,6,1);
% plot(visit_index,alpha_by_visit,'k-'); ylim([0 0.5]); box off; ylabel('\alpha'); xlabel('visits'); ylim([0 0.1]); box off;

figure(1); clf;

for mice = 1:25

p_reward    = zeros(6,numel(visit_index));
p_stay      = zeros(6,numel(visit_index));
yes_reward  = zeros(1,numel(visit_index));

p_reward(:,1)   = 1/6;
p_stay(:,1)     = 0;
alpha_by_visit  = (a_scale * exp(-(visit_index-1)./a_tau)) + 0.01;


    for vi = 2:numel(visit_index)
    
        % pass all estimates from previous step
        p_reward(:,vi)   = p_reward(:,vi-1);
        p_stay(:,vi)     = p_stay(:,vi-1);
    
        p_reward(p_reward(:,vi)<=0,vi)=0.001;
        p_reward(p_reward(:,vi)>1,vi)=1;
    
        % choose option
        if rand(1)<=epsilon
            checked_port(vi) = randperm(6,1);
        else
            checked_port(vi) = randsample(6,1,true,p_reward(:,vi));
        end
    
        % rewarded?
        if rand(1)<=reward_sched(checked_port(vi))
            yes_reward(vi)=1;
        end
    
        switch opto_on
            case 1
                if checked_port(vi)>2 && yes_reward(vi)==1
                    alpha_by_visit(vi)= alpha_by_visit(vi)*2;
                    title('Worst 4 bump');
                end
            case 2
                if yes_reward(vi)==1
                    alpha_by_visit(vi)= alpha_by_visit(vi)*2;
                    title('All bump');
                end
            case 3
                if checked_port(vi)>3 && yes_reward(vi)==0
                    alpha_by_visit(vi) = 0.01; % this is a way to make beta ~ 1
                    title('Unrewarded opto');
                end
        end
    
    
        switch decay_model
            case 'minus_alpha'
                beta(vi) = 1-alpha_by_visit(vi);
    
            case 'ind_static'
                beta(vi) = 0.975;
    
            case 'ind_dynamic'
                beta(vi) = 1-(2*alpha_by_visit(vi));
    
        end
    
        % update chosen option
        p_reward(checked_port(vi),vi)   = alpha_by_visit(vi)*yes_reward(vi) + beta(vi)*p_reward(checked_port(vi),vi-1);
        p_stay(checked_port(vi),vi)     = alpha_by_visit(vi)*yes_reward(vi) + beta(vi)*p_stay(checked_port(vi),vi-1);
    
        if live_plot || vi == numel(visit_index)
            % plot current estimates
            figure(1);
            subplot(1,6,1:2);
            for ports=1:6
                if mice==1
                    plot([1 vi],[reward_sched(ports) reward_sched(ports)],'-','Color',[port_colmap(ports,:)  0.1],'LineWidth',4); hold on;
                end
                plot(1:vi,p_reward(ports,1:vi),'-','Color',[port_colmap(ports,:) 0.2]);
            end
            axis([0 numel(visit_index) 0 1]); box off; xlabel('visits'); ylabel('P^{estim.} (rew)');
        end
    end
    
    for jj=1:6
        exp_rew_probs(jj) = sum(yes_reward(find(checked_port==jj))) ./ numel(find(checked_port==jj));
        chose_probs(jj) = numel(find(checked_port==jj)) ./ numel(visit_index);
    end
    
    % exp_rew_probs
    % chose_probs
    % reward_sched
    % p_reward(:,end)'
    
    figure(1);
    subplot(1,6,3:4);
    plot([0 1],[0 1],'k--'); hold on;
    quart_inds = [1:400;401:800;801:1200;1201:1600];
    for qi=4
        scatter(reward_sched,median(p_reward(:,quart_inds(qi,:)),2),25*qi,port_colmap,'filled','MarkerFaceAlpha',0.2);
    end
    axis([0 0.8 0 0.8]); box off; ylabel('P^{estim.} (rew)'); xlabel('P(rew)')
    
    figure(1);
    subplot(1,6,5:6);
    if mice==1
        hold off;
        plot([-5 0],[-5 0],'k--'); hold on;
    end
    quart_inds = [1:400;401:800;801:1200;1201:1600];
    for qi=4
        scatter(log(exp_rew_probs./sum(exp_rew_probs)),log(chose_probs./sum(chose_probs)),25*qi,port_colmap,'filled','MarkerFaceAlpha',0.2);
    end
    axis([-4 -0.5 -4 -0.5]); box off; xlabel('log reward'); ylabel('log choice');
end