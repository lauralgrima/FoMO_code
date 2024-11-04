% simple toy simulation of p_reward updates with both + / - arms
recalc_rew_seq = 1;
num_trials  = 400;
p_rew_alpha = zeros(1,num_trials);
p_rew_al_be = zeros(1,num_trials);
p_no_rew_al_be = zeros(1,num_trials);
p_rew_cum   = zeros(1,num_trials);
base_p      = zeros(1,num_trials);
% base_p(1:num_trials/3) = 0.75;
% base_p(num_trials/3:2*num_trials/2) = 0.5;
% base_p(2*num_trials/3:end) = 0.25;
base_p(1:num_trials/2) = 0.8;
base_p(num_trials/2:end) = 0.4;
alpha = ones(1,num_trials);
vis_inds = [1:200 1:200];
alpha=0.0033 + (0.1*exp(-vis_inds/75));
% figure(); plot(alpha)

sess_cmap = TNC_CreateRBColormap(8,'hue7');
beta  = 0.001;

if recalc_rew_seq
    rew         = zeros(1,num_trials);

    for jj=2:num_trials
        if rand(1)<=base_p(jj)
            rew(jj) = 1;
        else
            rew(jj) = -0.1;
        end
    end
end


for jj=2:num_trials

    % matching standard
    p_rew_cum(jj) = sum(rew(1:jj)./jj);

    % version that doesn’t need a buffer
    p_rew_alpha(jj) = alpha(end)*rew(jj) + (1-alpha(end))*p_rew_alpha(jj-1);

    if rew(jj)<0
        p_rew_al_be(jj) = (alpha(jj)*0.1)*rew(jj) + (1-alpha(jj))*p_rew_al_be(jj-1);
    else
        p_rew_al_be(jj) = alpha(jj)*rew(jj) + (1-alpha(jj))*p_rew_al_be(jj-1);
    end

end

    figure(1); hold off;
    plot(1:num_trials,base_p,'k-','color',[0.5 0.5 0.5],'LineWidth',4); hold on;
    plot(1:jj,p_rew_cum(1:jj),'k','LineWidth',2); hold on;
    plot(1:jj,p_rew_alpha(1:jj),'color',sess_cmap(3,:),'LineWidth',2); hold on;
    plot(1:jj,p_rew_al_be(1:jj),'color',sess_cmap(1,:),'LineWidth',2); hold on;
    % plot(1:jj,1+p_no_rew_al_be(1:jj),'g'); hold on;
    % plot(1:jj,((1+p_no_rew_al_be(1:jj))+p_rew_al_be(1:jj))./2,'r'); hold on;
    axis([0 num_trials 0 1]); box off;
    ylabel('P(rew) estimate');
    xlabel('Experience');
    drawnow; 
    legend({'True P(rew)','Matching','Static \alpha','Dynamic \alpha'})