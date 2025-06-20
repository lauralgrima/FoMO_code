visit_index = 1:1:1600;
a_tau = 200;
a_scale = 0.3;
alpha_by_visit = (a_scale * exp(-(visit_index-1)./a_tau)) + 0.005;

reward_sched = [0.59,0.5,0.37,0.36,0.17,0.15]

[port_colmap] = TNC_CreateRBColormap(6,'grima');

p_reward    = zeros(6,numel(visit_index));
p_stay      = zeros(6,numel(visit_index));
yes_reward  = zeros(1,numel(visit_index));

p_reward(:,1)   = 1/6;
p_stay(:,1)     = 0;

epsilon = 0.1;

opto_on = 1;

figure(1); clf;
subplot(1,6,1);
plot(visit_index,alpha_by_visit,'k-'); ylim([0 1]); box off;

for vi = 2:numel(visit_index)

    % pass all estimates from previous step
    p_reward(:,vi)   = p_reward(:,vi-1);
    p_stay(:,vi)     = p_stay(:,vi-1);

    p_reward(p_reward(:,vi)<0,vi)=0;

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
                alpha_by_visit(vi)= alpha_by_visit(vi)*1.5;
                title('Worst 4 bump');
            end
        case 2
            if yes_reward(vi)==1
                alpha_by_visit(vi)= alpha_by_visit(vi)*1.5;
                title('All bump');
            end
        case 3
            if checked_port(vi)>2 && yes_reward(vi)==0
                alpha_by_visit(vi)= alpha_by_visit(vi)*1.5;
                title('Unrewarded opto');
            end
    end

    % update chosen option
    p_reward(checked_port(vi),vi)   = alpha_by_visit(vi)*yes_reward(vi) + (1-alpha_by_visit(vi))*p_reward(checked_port(vi),vi-1);
    p_stay(checked_port(vi),vi)     = alpha_by_visit(vi)*yes_reward(vi) + (1-alpha_by_visit(vi))*p_stay(checked_port(vi),vi-1);

    % plot current estimates
    figure(1);
    subplot(1,6,2:5);
    hold off;
    for ports=1:6
        plot([1 vi],[reward_sched(ports) reward_sched(ports)],'-','Color',[port_colmap(ports,:)  0.05],'LineWidth',4); hold on;
        plot(1:vi,p_reward(ports,1:vi),'-','Color',[port_colmap(ports,:) 0.75]);
    end
    axis([0 numel(visit_index) 0 1]); box off;
end

figure(1);
subplot(1,6,6);
plot([0 1],[0 1],'k--'); hold on;
quart_inds = [1:400;401:800;801:1200;1201:1600];
for qi=1:4
    scatter(reward_sched,mean(p_reward(:,quart_inds(qi,:)),2),25*qi,port_colmap,'filled');
end
axis([0 1 0 1]); box off;