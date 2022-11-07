
%% Exploring the licking data a bit

lk_kern = TNC_CreateGaussian(500,200,1000,1);
lk_kern = lk_kern./ max(lk_kern);
% lk_kern(1:500)=0;
% lk_kern = [0 ones(1,1000) 0];

% make a matrix of dimension 6 ports x total time
lk_6p6_raw = csvread('lick_df.csv',1,1);
hexa.all_raw = zeros(6,lk_6p6_raw(end,2));
lk_6p6_vid = csvread('video_data_df.csv',1,1);

prior_port_ids = [lk_6p6_raw(1,3) ; lk_6p6_raw(1:end-1,3)]';
for qq=1:6
    hexa.all_raw(qq,lk_6p6_raw(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1,1)) = 1;
    % find the prior port
    hexa.prior_port(qq,lk_6p6_raw(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1,1)) = prior_port_ids(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1);
end

figure(11); clf;
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
end

eg_port = 2;
figure(12); clf;
subplot(2,4,1:2);
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
end
box off;

lk_intervals = [0 diff(find(hexa.all_raw(eg_port,:)==1))];
lk_times = find(hexa.all_raw(eg_port,:)==1);

subplot(2,4,5:6);
plot(conv(hexa.all_raw(eg_port,:),lk_kern,'same'));
plot(lk_times(lk_intervals>1000),ones(1,sum(lk_intervals>1000)),'*');
box off;

subplot(2,4,[3 7]);
hist(log10(diff(find(hexa.all_raw(eg_port,:)==1))),250);


% add a piece of logic to check for switched port
port_ids = hexa.all_raw .* [1:6]';

figure(13); clf;
for qq=1:6
    lk_intervals = [0 diff(find(hexa.all_raw(qq,:)==1))];
    port_changes = [hexa.prior_port(qq,find(hexa.all_raw(qq,:)==1))];
    visit_cnts(qq,1) = sum(lk_intervals>1000);
    visit_cnts(qq,2) = sum(port_changes~=qq);
    subplot(1,6,qq);
    swarmchart(port_changes,log10(lk_intervals),'.'); hold on;
    plot([0 7],[3 3],'k-');
    axis([0 7 0 7]);
    xlabel('Previous port'); ylabel('Lick interval (log ms)'); title(['port: ' num2str(qq)]);
end
figure(12);
subplot(2,4,[4 8]);
plot(1:6,sum(hexa.all_raw,2),'o-');
yyaxis right
plot(1:6,visit_cnts(:,1),'o-');
plot(1:6,visit_cnts(:,2),'o-');


%% Analyzing visits at video rate
% 
lk_6p6_raw = csvread('lick_df.csv',1,1);
lk_6p6_vid = csvread('video_data_df.csv',1,1);

lk_kern = TNC_CreateGaussian(500,200,1000,1);
lk_kern = lk_kern./ max(lk_kern);

[pos] = TNC_BigBoiConverter(lk_6p6_vid);

hexa.all_raw = zeros(6,max(lk_6p6_vid(:,1)));
hexa.prior_port = zeros(6,max(lk_6p6_vid(:,1)));
hexa.pos(1,pos.frms) = pos.xf;
hexa.pos(2,pos.frms) = pos.yf;
hexa.vel = medfilt1(pos.v);

prior_port_ids = [lk_6p6_raw(1,3) ; lk_6p6_raw(1:end-1,3)]';
for qq=1:6

    % find samples for port==qq
    samps = find(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1);
    frames = lk_6p6_raw(samps,5);
    hexa.all_raw(qq,frames(frames>0)) = 1;
    % find the prior port
    hexa.prior_port(qq,frames(frames>0)) = prior_port_ids(samps(frames>0));
end

figure(11); clf;
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
end

eg_port = 2;
figure(12); clf;
subplot(2,4,1:2);
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
end
box off;

lk_intervals = [0 diff(find(hexa.all_raw(eg_port,:)==1))];
lk_times = find(hexa.all_raw(eg_port,:)==1);

subplot(2,4,5:6);
plot(conv(hexa.all_raw(eg_port,:),lk_kern,'same'));
plot(lk_times(lk_intervals>200),ones(1,sum(lk_intervals>200)),'*');
box off;

subplot(2,4,[3 7]);
hist(log10(diff(find(hexa.all_raw(eg_port,:)==1))),250);

% add a piece of logic to check for switched port
port_ids = hexa.all_raw .* [1:6]';

visit_iti_thresh = 400;
hexa.visits = zeros(size(hexa.all_raw));

% Not quite finished, but need to know which vists were rewarded
% hexa.rewards = zeros(1,size(hexa.all_raw,2));
% hexa.rewards = 

figure(13); clf;
for qq=1:6
    % time between licks on same port
    lk_intervals = [0 diff(find(hexa.all_raw(qq,:)==1))];
    % port change counts as new visit by definition
    port_changes = [hexa.prior_port(qq,find(hexa.all_raw(qq,:)==1))];
    % compute distance traveled between licks
    all_licks = find(hexa.all_raw(qq,:)==1);
    dist_iti = zeros(size(all_licks));
    dist_iti(1) = 1000;
    for mm=2:numel(all_licks)
        dist_iti(mm) = trapz(hexa.vel(all_licks(mm-1):all_licks(mm)));
    end

    visit_cnts(qq,1) = sum(lk_intervals>visit_iti_thresh);
    visit_cnts(qq,2) = sum(port_changes~=qq);

    visit_true = find(lk_intervals>visit_iti_thresh | dist_iti>500 | port_changes~=qq);
    all_licks = find(hexa.all_raw(qq,:)==1);
    hexa.visits(qq,all_licks(visit_true)) = 1;
    
    subplot(1,6,qq);
    swarmchart(port_changes,log10(lk_intervals),dist_iti,'.'); hold on;
    plot([0 7],[log10(visit_iti_thresh) log10(visit_iti_thresh)],'k-');
    axis([0 7 0 7]);
    xlabel('Previous port'); ylabel('Lick interval (log ms)'); title(['port: ' num2str(qq)]);
end

figure(12);
subplot(2,4,[4 8]);
plot(1:6,sum(hexa.all_raw,2),'o-'); hold on;
yyaxis right
plot(1:6,visit_cnts(:,1),'o-');
plot(1:6,visit_cnts(:,2),'o-');

figure(21);
subplot(1,5,1:4);
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
    plot(hexa.visits(qq,:)+qq,'k','linewidth',1); hold on;
end
box off;
subplot(1,5,5);
hexa.visits_iti = diff(find(sum(hexa.visits,1)==1));
hist(log10(hexa.visits_iti),50);
box off;
xlabel('Visit Intervals (log ms)');
ylabel('Count');

figure(20); clf;
plot(hexa.pos(1,:),hexa.pos(2,:),'color',[0 0 0 0.1]);
pmap = TNC_CreateRBColormap(6,'mapb');
for eg_port = 1:6
    lk_times = find(hexa.all_raw(eg_port,:)==1);
    hold on;
    plot(hexa.pos(1,lk_times),hexa.pos(2,lk_times),'.','color',pmap(eg_port,:));
end

%% Examine empirical hazard based p_reward estimates from data
max_tsteps = size(hexa.visits,2);
p_reward = zeros(size(hexa.visits));
p_NOreward = zeros(size(hexa.visits));
ex_map = TNC_CreateRBColormap(100,'exag');
cat_map = TNC_CreateRBColormap(8,'mapb');

all_rew_inds = find(lk_6p6_raw(:,4)==1);
all_rew_ports = lk_6p6_raw(all_rew_inds,3);
all_rew_frames = lk_6p6_raw(all_rew_inds,5);

for zzz=1:6
    hexa.rew_port_v(zzz).fr = all_rew_frames(all_rew_ports==zzz)';
end

figure(20); clf;
for zzz=1:6
    for t=1:max_tsteps
                   
        all_rew_t = hexa.rew_port_v(zzz).fr(hexa.rew_port_v(zzz).fr<=t & hexa.rew_port_v(zzz).fr>0);
        all_vis_t = find(hexa.visits(zzz,1:t)==1);

        if numel(all_rew_t)>0
        
            if numel(find(t==all_rew_t))==1
                pdf                = hist([all_rew_t(1) diff(all_rew_t)],0:1:1e6);
                cdf                = cumsum(pdf)./sum(pdf);
%                 cdf                = cumsum(pdf)./numel(find(hexa.visits(zzz,1:t)==1));

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

    plot(p_reward(zzz,:)+zzz,'color',cat_map(zzz,:)); 
    hold on; plot(all_rew_t-1,p_reward(zzz,all_rew_t-1)+zzz,'*','color',cat_map(zzz,:)); 
    hold on; plot(find(hexa.visits(zzz,:)==1)-1,-p_NOreward(zzz,find(hexa.visits(zzz,:)==1)-1)+zzz,'.','color',cat_map(zzz,:)); 
    plot(-p_NOreward(zzz,:)+zzz,'color',cat_map(zzz,:)./2); 
    axis([0 max_tsteps 1 7]); box off;
    drawnow;

end


%% Examine p_hazard for each port choice
figure(21);
for zzz=1:6
    all_rew_t = hexa.rew_port_v(zzz).fr(hexa.rew_port_v(zzz).fr<t & hexa.rew_port_v(zzz).fr>0);
    [sink] = TNC_ExtTrigWins3d(p_reward,find(hexa.visits(zzz,:)==1)-1,[0 0]);
%     [sink] = TNC_ExtTrigWins3d(p_reward,all_rew_t-1,[0 0]);
    subplot(1,6,zzz);
    imagesc(sink.win_avg,[0 max(sink.win_avg)]); colormap(ex_map); title(zzz);
    axis off;
end


%%
figure(22); clf;
xvals = 0:1:5e5;
subplot(121);
    for zzz=1:6
        all_rew_t = hexa.rew_port_v(zzz).fr(hexa.rew_port_v(zzz).fr<t & hexa.rew_port_v(zzz).fr>0);    
        tmp = hist([all_rew_t(1) diff(all_rew_t)],xvals);
        cdf = cumsum(tmp);
        thresh = find(cdf./cdf(end)>1/6,1);
        xvals(thresh);        
        semilogx(cdf./cdf(end)); hold on;
        set(gca,'XLim',[50 5e5]);
    end
    
    plot([1 10e5],[1 1]/6,'k--');
    title('Reward intervals'); box off;

subplot(122);
    for zzz=1:6
        all_rew_t = find(hexa.visits(zzz,:)==1);   
        tmp = hist(diff(all_rew_t),xvals);
        cdf = cumsum(tmp);
%         cdf = cdf./cdf(end);
        semilogx(cdf./cdf(end)); hold on;
        set(gca,'XLim',[50 5e5]);
    end
    title('Visit intervals'); box off;

%% Get experienced reward intervals on each port

all_rew_inds = find(lk_6p6_raw(:,4)==1);
all_rew_ports = lk_6p6_raw(all_rew_inds,3);
all_rew_frames = lk_6p6_raw(all_rew_inds,1);
figure(57); clf;
cat_map = TNC_CreateRBColormap(8,'mapb');

for qq=1:6

    hexa.port_rew(qq).ts = all_rew_frames(all_rew_ports==qq);
    hexa.port_rew(qq).iti = diff(hexa.port_rew(qq).ts)/1000;
    disp(['Port ' num2str(qq) ' mean iti: ' num2str(mean(hexa.port_rew(qq).iti)) ' +/- ' num2str(std(hexa.port_rew(qq).iti))]);
    hexa.port_rew(qq).histX = 0.01:0.1:4;
    hexa.port_rew(qq).histY = hist(log10(hexa.port_rew(qq).iti),hexa.port_rew(qq).histX);
    plot(hexa.port_rew(qq).histX , hexa.port_rew(qq).histY , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;

end

xlabel('Log10 Interreward Interval'); ylabel('Count'); legend; box off;

%% Goal is to explain the hexa.visits data matrix

hexa.port_chk_t = find(sum(hexa.visits,1)==1);
hexa.port_chk   = hexa.visits(:,hexa.port_chk_t);

tmap = TNC_CreateRBColormap(8,'exag');
hexa.trans_mat = zeros(6,6);

% compute the transition matrix of visits
for qq=1:6
    tmp = find(hexa.port_chk(qq,1:end-1)==1);
    trans_cnt = zeros(1,6);
    for pp=tmp
        next_visit = find(hexa.port_chk(:,pp+1)==1);
        trans_cnt(next_visit) = trans_cnt(next_visit) + 1;
    end
    hexa.trans_mat(qq,:) = trans_cnt;
%     hexa.trans_mat(:,qq) = trans_cnt';
end

% Want to compare these two matrices + reward fraction with models 
figure(31); clf;
subplot(1,4,1:2);
imagesc(hexa.port_chk);
xlabel('Unique visits'); ylabel('Port');
subplot(1,4,3);
plot(sum(hexa.port_chk,2),1:6,'o-');
xlabel('Total visits'); ylabel('Port');
axis([0 size(hexa.port_chk,2)/2 1 6]); box off;
set(gca,'YDir','reverse');

subplot(1,4,4);
imagesc(hexa.trans_mat./size(hexa.port_chk,2),[0 0.15]); colormap(tmap); colorbar;
xlabel('Port'); ylabel('Port');
title('Transition probability');

%% Model comparisons test script 

clear hexa_model;
hexa_model.seed = randperm(100,1);
rng(hexa_model.seed);

hexa_model.interportdist = ...
[0	14	18	70	72.2	65.5 ;...
14	0	22.8	56	65.5	42;...
18	22.8	0	72.2	70	56;...
70	56	72.2	0	18	22.8;...
72.2	65.5	70	18	0	14;...
65.5	42	56	22.8	14	0];

% sampling rate is now set to camera frame rate
frame_rate = 200;

% Make the BigBoi reward environment
hexa_model.rew_sched = zeros(size(hexa.visits));
hexa_model.rew_sched(1,1:30*frame_rate:end) = 1;
hexa_model.rew_sched(2,200:60*frame_rate:end) = 1;
hexa_model.rew_sched(3,400:240*frame_rate:end) = 1;
hexa_model.rew_sched(4,600:240*frame_rate:end) = 1;
hexa_model.rew_sched(5,800:1200*frame_rate:end) = 1;
hexa_model.rew_sched(6,1000:2400*frame_rate:end) = 1;
hexa_model.rew_sched(:,1) = 1;

max_reward = sum(sum(hexa_model.rew_sched));

% Pass in data file and policy choice
policy.type = 'e-proportional'; % out of type = {'softmax','greedy','e-greedy','random','proportional','e-proportional'}
policy.params.epsilon = 0.075;

belief.type = 'matchP-shift-spatial'; % out of type = {'win-stay','proportional','kernel','spatial','pdf','pdf-space'}
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

max_tsteps = size(hexa.visits,2);
sample_logic = sum(hexa.visits,1);
hexa_model.visits = zeros(size(hexa.visits));
hexa_model.rewards = zeros(size(hexa.visits));
reward_available = zeros(size(hexa.visits));
p_reward = zeros(size(hexa.visits));
p_reward(:,1) = 1/6;

for t=2:max_tsteps-1
    
   reward_available(reward_available(:,t)==0,t) = hexa_model.rew_sched(reward_available(:,t)==0,t);
   reward_available(:,t+1) = reward_available(:,t);

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

           case 'hazard' %- attempt to estimate true posterior P(rew|port,t)
               % at this moment in time what is P(rew) over ports based
               % upon an estimate of probability over time (i.e. like the
               % hazard)
               for zzz=1:6
                   
                   all_rew_t = find(hexa_model.rewards(zzz,1:t)==1);
                   
                   if numel(all_rew_t)>1
                       
                       hexa_model.last_rew(zzz) = all_rew_t(end);
                       hexa_model.mu_rew(zzz) = mean([all_rew_t(1) diff(all_rew_t)]);
                       hexa_model.sg_rew(zzz) = std([all_rew_t(1) diff(all_rew_t)]);
                       if hexa_model.sg_rew(zzz)<hexa_model.mu_rew(zzz)*0.2
                           hexa_model.sg_rew(zzz)=hexa_model.mu_rew(zzz)*0.2; % standard scalar timing
                       end

                       % compute pdf
                       pdf = TNC_CreateGaussian(hexa_model.mu_rew(zzz),hexa_model.sg_rew(zzz),t-hexa_model.last_rew(zzz),1);

                       % integrate pdf from t_prev_rew:t to get p_reward(t)
                       p_reward(zzz,t) = trapz(pdf);

                   else
                       
                       hexa_model.last_rew(zzz) = NaN;
                       hexa_model.mu_rew(zzz) = NaN;
                       hexa_model.sg_rew(zzz) = NaN;
                       p_reward(zzz,t) = policy.params.epsilon;

                   end
               end


           otherwise % do nothing

       end
    
       % Combine beliefs into a meta belief (discover apprpriate weightings)


   end
    
end

disp(['Model rewards collected: ' num2str(sum(sum(hexa_model.rewards))) ' ; ' num2str(100*(sum(sum(hexa_model.rewards)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Mouse rewards collected: ' num2str(sum(lk_6p6_raw(:,4))) ' ; ' num2str(100*(sum(lk_6p6_raw(:,4)))/(sum(sum(hexa_model.rew_sched)))) '%'])
disp(['Max rewards available: ' num2str(sum(sum(hexa_model.rew_sched)))])

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
plot(sum(hexa.port_chk,2),1:6,'ko-');
legend({'Model','Data'});
xlabel('Total visits'); ylabel('Port');
set(gca,'YDir','reverse');
axis([0 size(hexa_model.port_chk,2)/2 1 6]); box off;
subplot(1,4,4);
imagesc(hexa_model.trans_mat./size(hexa_model.port_chk,2),[0 0.15]); colormap(tmap); colorbar;
xlabel('Port'); ylabel('Port');
title('Transition probability');

