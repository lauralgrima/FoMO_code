% simple toy simulation of p_reward updates with both + / - arms

num_trials  = 1500;
p_rew_alpha = zeros(1,num_trials);
p_rew_al_be = zeros(1,num_trials);
p_no_rew_al_be = zeros(1,num_trials);
p_rew_cum   = zeros(1,num_trials);
rew         = zeros(1,num_trials);
base_p      = zeros(1,num_trials);
base_p(1:num_trials/3) = 0.67;
base_p(num_trials/3:2*num_trials/2) = 0.63;
base_p(2*num_trials/3:end) = 0.15;
alpha = 0.1;
beta  = 0.01;

for jj=2:num_trials

    if rand(1)<=base_p(jj)
        rew(jj) = 1;
    end

    % matching standard
    p_rew_cum(jj) = sum(rew(1:jj)./jj);

    % version that doesn’t need a buffer
    p_rew_alpha(jj) = alpha*rew(jj) + (1-alpha)*p_rew_alpha(jj-1);

    p_rew_al_be(jj) = alpha*rew(jj) + (1-alpha)*p_rew_al_be(jj-1);

    % version that learns from rew and 
    % beta = 0.05*(1-p_rew_al_be(jj-1));
    p_no_rew_al_be(jj) = beta*(rew(jj)-1) + (1-beta)*p_no_rew_al_be(jj-1);

    figure(1); hold off;
    plot(1:num_trials,base_p,'k--'); hold on;
    plot(1:jj,p_rew_cum(1:jj),'k'); hold on;
    % plot(1:jj,p_rew_alpha(1:jj),'r'); hold on;
    plot(1:jj,p_rew_al_be(1:jj),'b'); hold on;
    plot(1:jj,1+p_no_rew_al_be(1:jj),'g'); hold on;
    plot(1:jj,((1+p_no_rew_al_be(1:jj))+p_rew_al_be(1:jj))./2,'r'); hold on;
    axis([0 num_trials 0 1]); box off;
    drawnow; 
end