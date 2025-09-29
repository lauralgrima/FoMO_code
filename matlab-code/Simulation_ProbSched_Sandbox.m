
visit_index = 1:1:2000;
a_tau = 100;
a_scale = 0.0099;
alpha_by_visit = (a_scale * exp(-(visit_index-1)./a_tau)) + 0.0033;
alpha_by_visit_stim = (a_scale * exp(-(visit_index-1)./800)) + 0.0033;

clear Summary_UpdateRule_OptoSim

% figure(1); clf; subplot(121); plot(alpha_by_visit); hold on; plot(alpha_by_visit_stim)

reward_sched = [0.59,0.5,0.37,0.36,0.17,0.15];

[port_colmap] = TNC_CreateRBColormap(6,'grima');

sensitivity = [];

p_reward    = zeros(6,numel(visit_index));
p_stay      = zeros(6,numel(visit_index));
yes_reward  = zeros(1,numel(visit_index));

p_reward(:,1)   = 1/6;
p_stay(:,1)     = 0;

epsilon = 0.1;

live_plot = 0;
decay_model = 'minus_alpha'

for opto_on = [0 -1 2 4]

    figure(7+opto_on); clf;
    
    for mice = 1:25
    
        p_reward    = zeros(6,numel(visit_index));
        p_stay      = zeros(6,numel(visit_index));
        yes_reward  = zeros(1,numel(visit_index));
        
        switch opto_on
            case 2
                p_reward(:,1)   = flip(reward_sched);
                p_stay(:,1)     = 0;
            case 4
                p_reward(:,1)   = flip(reward_sched);
                p_stay(:,1)     = 0;
            otherwise
                p_reward(:,1)   = 1/6;
                p_stay(:,1)     = 0;
        end

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
    
            this_alpha = alpha_by_visit(vi);
            
            switch opto_on
                case -1
                    if checked_port(vi)>2 && yes_reward(vi)==1
                        this_alpha = alpha_by_visit_stim(vi);
                        title('Worst 4 bump (tau)');
                    end
                case 0
                    title('Control simulation');            
                    
                case 1
                    if checked_port(vi)>2 && yes_reward(vi)==1
                        this_alpha = alpha_by_visit(vi)*1.5;
                        title('Worst 4 bump');
                    end
                case 2
                    if yes_reward(vi)==1
                        this_alpha = alpha_by_visit_stim(vi);
                        title('All bump');
                    end
                case 3
                    if checked_port(vi)>3 && yes_reward(vi)==0
                        this_alpha = 0.01; % this is a way to make beta ~ 1
                        title('Unrewarded opto');
                    end
                case 4
                    title('Control switch simulation');            
            end
        
        
            switch decay_model
                case 'minus_alpha'
                    this_beta = 1-this_alpha;
        
                case 'ind_static'
                    this_beta = 0.975;
        
                case 'ind_dynamic'
                    this_beta = 1-(2*this_alpha);
        
            end
        
            % update chosen option
            p_reward(checked_port(vi),vi)   = this_alpha*yes_reward(vi) + this_beta*p_reward(checked_port(vi),vi-1);
            p_stay(checked_port(vi),vi)     = this_alpha*yes_reward(vi) + this_beta*p_stay(checked_port(vi),vi-1);
        
            if live_plot || vi == numel(visit_index)
                % plot current estimates
                figure(7+opto_on);
                subplot(1,6,1:2);
                for ports=1:6
                    if mice==1
                        plot([numel(visit_index) numel(visit_index)+50],[reward_sched(ports) reward_sched(ports)],'-','Color',[port_colmap(ports,:)],'LineWidth',4); hold on;
                    end
                    plot(1:vi,p_reward(ports,1:vi),'-','Color',[port_colmap(ports,:)]);
                end
                axis([0 numel(visit_index)+50 0 1]); box off; xlabel('visits'); ylabel('P^{estim.} (rew)');
            end
        end
        
        for jj=1:6
            valid_reward    = yes_reward(round(end/2):end);
            valid_port      = checked_port(round(end/2):end);
            exp_rew_probs(jj) = sum(valid_reward(find(valid_port==jj))) ./ numel(find(valid_port==jj));
            chose_probs(jj) = numel(find(valid_port==jj)) ./ numel(valid_port);
        end
    
        sens_fit = polyfit(log(exp_rew_probs./sum(exp_rew_probs)),log(chose_probs./sum(chose_probs)),1);
        sensitivity(mice) = sens_fit(1);
        
        figure(7+opto_on);
        subplot(1,6,3:4);
        plot([0 1],[0 1],'k--'); hold on;
        quart_inds = [1:400;401:800;801:1200;1201:1600];
        for qi=4
            scatter(reward_sched,median(p_reward(:,quart_inds(qi,:)),2),25*qi,port_colmap,'filled');
        end
        axis([0 0.8 0 0.8]); box off; ylabel('P^{estim.} (rew)'); xlabel('P(rew)')
        
        figure(7+opto_on);
        subplot(1,6,5:6);
        if mice==1
            hold off;
            plot([-5 0],[-5 0],'k--'); hold on;
        end
        quart_inds = [1:400;401:800;801:1200;1201:1600];
        for qi=4
            scatter(log10(exp_rew_probs./sum(exp_rew_probs)),log10(chose_probs./sum(chose_probs)),25*qi,port_colmap,'filled');
        end
        axis([-1.5 0 -1.5 0]); box off; xlabel('log reward'); ylabel('log choice');
    
        Summary_UpdateRule_OptoSim.reward_sched             = reward_sched;
        % Summary_UpdateRule_OptoSim.log_reward(mice,:)       = log10(exp_rew_probs./sum(exp_rew_probs));
        % Summary_UpdateRule_OptoSim.log_choice(mice,:)       = log10(chose_probs./sum(chose_probs));
        % Summary_UpdateRule_OptoSim.sims(mice).reward        = yes_reward;
        % Summary_UpdateRule_OptoSim.sims(mice).p_reward      = p_reward;
        % Summary_UpdateRule_OptoSim.sims(mice).checked_port  = checked_port;

        switch opto_on
            
            case -1
                Summary_UpdateRule_OptoSim.rew_prob_observed_4worst(mice,:) = exp_rew_probs;
                Summary_UpdateRule_OptoSim.chosen_port_probs_4worst(mice,:) = chose_probs;

            case 2
                Summary_UpdateRule_OptoSim.rew_prob_observed_all6(mice,:) = exp_rew_probs;
                Summary_UpdateRule_OptoSim.chosen_port_probs_all6(mice,:) = chose_probs;

            case 0
                Summary_UpdateRule_OptoSim.rew_prob_observed(mice,:) = exp_rew_probs;
                Summary_UpdateRule_OptoSim.chosen_port_probs(mice,:) = chose_probs;

        end
    
    end

    figure(7+opto_on);
    subplot(1,6,5:6);
    text(-1.5,-0.2,['s = ' num2str(mean(sensitivity)) ' +/- ' num2str(std(sensitivity))]);

    final_boxplot(:,find(opto_on==[0 -1 4 2])) = sensitivity;

end

figure(1); clf; 
subplot(211); boxchart(final_boxplot(:,1:2));
xticklabels({'Control' 'Stim 3-6'}); ylim([0.2 1]); ylabel('sensitivity');
subplot(212); boxchart(final_boxplot(:,3:4));
xticklabels({'Switch Control' 'Switch Stim All'}); ylim([0 0.6]); ylabel('sensitivity');
