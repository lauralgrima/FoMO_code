
%% Filenames
path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';

% filenames = {'6PG5_NAc_conc_beh.csv'};
% filenames = {'6PG9_DMS_conc_beh.csv'};
filenames = {'6PG12_NAc_conc_beh.csv'};

%% Core workflow:

model_dist = 1; reps = 20; clear hexa_model_an;
for mm = 1
    for session = 1
        [hexa_data]     = HX_load_csv([path filenames{mm}], 0, 1);
        [hexa_data_an]  = HX_analyze_session(hexa_data,session,1);
        if model_dist
            pfm = figure(61); clf;
            for jj=1:reps
                [hexa_model]    = HX_model_session(hexa_data_an,'e-proportional','matching',0,0);
                hexa_model_an.sim_reps.ideal(jj,:)      = cumsum(sum(hexa_model.ideal,1));
                hexa_model_an.sim_reps.random(jj,:)     = cumsum(sum(hexa_model.random,1));
                hexa_model_an.sim_reps.rewards(jj,:)    = cumsum(sum(hexa_model.rewards,1));                
                hexa_model_an.sim_reps.Sideal(jj,:)     = hexa_model.slope.ideal;
                hexa_model_an.sim_reps.Srandom(jj,:)    = hexa_model.slope.random;
                hexa_model_an.sim_reps.Srewards(jj,:)   = hexa_model.slope.model;    
                hexa_model_an.sim_reps.visits(jj,:)     = sum(hexa_model.visits,2);
            end
            subplot(131); hold off;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.ideal,1),2.*std(hexa_model_an.sim_reps.ideal,1),{'color',[0 0 0]}); hold on;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.random,1),2.*std(hexa_model_an.sim_reps.random,1),{'color',[0.5 0.5 0.5]}); hold on;
            shadedErrorBar(1:size(hexa_model_an.sim_reps.ideal,2),mean(hexa_model_an.sim_reps.rewards,1),2.*std(hexa_model_an.sim_reps.rewards,1),{'color',[0.5 0 0.16]}); hold on;
            plot(1:size(hexa_model_an.sim_reps.ideal,2),cumsum(sum(hexa_data_an.rewards,1)),'-','color',[1 0 0.33],'linewidth',2); hold on;
            box off;
            ylabel('Cumulative rewards'); xlabel('Unique Port Visits');
            axis([0 size(hexa_model.ideal,2) 0 max(cumsum(sum(hexa_model.ideal,1)))]);

            subplot(132); hold off;
            shadedErrorBar(hexa_model.slope.x,mean(hexa_model_an.sim_reps.Sideal,1),2.*std(hexa_model_an.sim_reps.Sideal,1),{'color',[0 0 0]}); hold on;
            shadedErrorBar(hexa_model.slope.x,mean(hexa_model_an.sim_reps.Srandom,1),2.*std(hexa_model_an.sim_reps.Srandom,1),{'color',[0.5 0.5 0.5]}); hold on;
            shadedErrorBar(hexa_model.slope.x,mean(hexa_model_an.sim_reps.Srewards,1),2.*std(hexa_model_an.sim_reps.Srewards,1),{'color',[0.5 0 0.16]}); hold on;
            plot(hexa_model.slope.x,hexa_model.slope.mouse,'-','color',[1 0 0.33],'linewidth',2); hold on;
            box off;
            ylabel('Local Income'); xlabel('Session Time');
            axis([0 size(hexa_model.ideal,2) 0.015 0.075]);

            subplot(133); hold off;
            shadedErrorBar(1:6,mean(hexa_model_an.sim_reps.visits,1),2.*std(hexa_model_an.sim_reps.visits,1),{'color',[0.5 0 0.16]}); hold on;
            plot(1:6,sum(hexa_data_an.visits,2),'o-','color',[1 0 0.33],'linewidth',2);
            xlabel('Total visits'); ylabel('Port');
            axis([0 7 0 max(sum(hexa_data_an.visits,2))*1.5]); 
            box off;

        else
            [hexa_model]    = HX_model_session(hexa_data_an,'e-proportional','match-shift',1,1);
        end
    end
end

%% plot HX_model_session ITI histograms

cat_map = TNC_CreateRBColormap(8,'mapb');
figure(107); clf;
[all_rew_ports,all_rew_inds] = find(hexa_model.rewards==1);
[all_vis_ports,all_vis_inds] = find(hexa_model.visits==1);

for qq=1:6
    hexa_model_an.port_rew(qq).ts = all_rew_inds(all_rew_ports==qq);
    hexa_model_an.port_rew(qq).iti = diff(hexa_model_an.port_rew(qq).ts);
    disp(['Port ' num2str(qq) ' mean ri: ' num2str(mean(hexa_model_an.port_rew(qq).iti)) ' +/- ' num2str(std(hexa_model_an.port_rew(qq).iti))]);
    hexa_model_an.port_rew(qq).histX = 0.1:0.1:4;
    hexa_model_an.port_rew(qq).histY = hist(log10(hexa_model_an.port_rew(qq).iti),hexa_model_an.port_rew(qq).histX);
    subplot(131);
    plot(hexa_model_an.port_rew(qq).histX , cumsum(hexa_model_an.port_rew(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    plot(hexa_data_an.port_rew(qq).histX , cumsum(hexa_data_an.port_rew(qq).histY) , '-.', 'color' , [cat_map(qq,:) 0.5] , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','1d','2','2d','3','3d','4','4d','5','5d','6','6d'},'Location','northwest');
        xlabel('Log10 Interreward Interval'); ylabel('Count'); box off;
    end

    subplot(132)
    hexa_model_an.port_vis(qq).ts = all_vis_inds(all_vis_ports==qq);
    hexa_model_an.port_vis(qq).iti = diff(hexa_model_an.port_vis(qq).ts);
    disp(['Port ' num2str(qq) ' mean vi: ' num2str(mean(hexa_model_an.port_vis(qq).iti)) ' +/- ' num2str(std(hexa_model_an.port_vis(qq).iti))]);
    hexa_model_an.port_vis(qq).histX = 0.1:0.1:4;
    hexa_model_an.port_vis(qq).histY = hist(log10(hexa_model_an.port_vis(qq).iti),hexa_model_an.port_vis(qq).histX);

    plot(hexa_model_an.port_vis(qq).histX , cumsum(hexa_model_an.port_vis(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    plot(hexa_data_an.port_vis(qq).histX , cumsum(hexa_data_an.port_vis(qq).histY) , '-.', 'color' , [cat_map(qq,:) 0.5], 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        % legend({'1','1d','2','2d','3','3d','4','4d','5','5d','6','6d'},'Location','northwest');
        xlabel('Log10 Intervisit Interval'); ylabel('Count'); box off;
    end
    
    
end

%% Calculate a transition matrix from visits matrix

Nback = 1;

tmp = find(sum(hexa_model.visits,1)==1);
[~,visit_list] = max(hexa_model.visits(:,tmp),[],1);
% [~,visit_list] = max(hexa_model.visits(:,tmp(1:round(end/3))),[],1);
% [~,visit_list] = max(hexa_model.visits(:,tmp(round(end/3):end)),[],1);


[trans_mat] = HX_ComputeTransitionMatrix(visit_list,25,Nback);
title(['SIMULATION; Nback=' num2str(Nback) ' trans. matrix']);


tmp = find(sum(hexa_data_an.visits,1)==1);
[~,visit_list_data] = max(hexa_data_an.visits(:,tmp),[],1);
% [~,visit_list_data] = max(hexa_data_an.visits(:,tmp(1:round(end/3))),[],1);
% [~,visit_list_data] = max(hexa_data_an.visits(:,tmp(round(end/3):end)),[],1);

[trans_mat_data] = HX_ComputeTransitionMatrix(visit_list_data,26,Nback);
title(['DATA; Nback=' num2str(Nback) ' trans. matrix']);


%% Data analysis questions:

% does the p(return|reward) evolve with a different timecourse to
% p(visit|reward history) across animals?

% does the time course of p(return|reward) correlate with / explain income
% slope inflections?

filenames   = dir();
Nback       = 1;

for mm = 1
    for session = 1

        [hexa_data]     = HX_load_csv([path filenames{mm}], 0, 1);
        [hexa_data_an]  = HX_analyze_session(hexa_data,session,1);

        [trans_mat_data] = HX_ComputeTransitionMatrix(visit_list_data,26,Nback);
        title(['DATA; Nback=' num2str(Nback) ' trans. matrix']);    
        
    end
end

%% Examine new 6 state markov version

% filenames = {'6PG5_NAc_conc_beh.csv'};
filenames = {'6PG12_NAc_conc_beh.csv'};

path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';
mm = 1; session = 1;
[hexa_data]     = HX_load_csv([path filenames{mm}], 0, 1);
[hexa_data_an]  = HX_analyze_session(hexa_data,session,1);

[hexa_model]    = HX_model_session_2(hexa_data_an,'e-proportional','p_check_match',-1,0);

cat_map = TNC_CreateRBColormap(8,'cat2');
cat_map = cat_map([1 3:6 8],:);
trial_win = 2;
trial_kernel = TNC_CreateGaussian(500,trial_win,1000,1);

figure(200); clf;
scatter(hexa_data_an.da_resp_all.t,hexa_data_an.da_resp_all.r,50,hexa_data_an.da_resp_all.p,'filled','MarkerFaceAlpha',0.25); colormap(cat_map); hold on;
plot(hexa_data_an.da_resp_all.t,conv( hexa_data_an.da_resp_all.r  , trial_kernel , 'same' ) ,'color' , [0 0.67 1 0.5], 'linewidth' , 4);
axis([0 max(visit_indices) -500 1500]);
yyaxis right;
visit_indices       = find(hexa_data.unique_vis==1 & ~isnan(hexa_data.photo_i) & hexa_data.session_n==session);
rew_visit_ids       = hexa_data.rewarded(visit_indices);
plot(hexa_data.event_time(visit_indices(rew_visit_ids==1)),conv( [1 diff(hexa_data.event_time(visit_indices(rew_visit_ids==1)))']  , trial_kernel , 'same' ) ,'color' , [0 0 0 0.5], 'linewidth' , 4);

figure(201); clf;
clear tmp_cross_corr
shifts = 50;
rew_rate = conv( [1 diff(hexa_data.event_time(visit_indices(rew_visit_ids==1)))']  , trial_kernel , 'same' );
da_resp_mag = conv( hexa_data_an.da_resp_all.r  , trial_kernel , 'same' );
cross_corr = xcorr(rew_rate-mean(rew_rate),da_resp_mag-mean(da_resp_mag),shifts );
num_mc = 200;
for kk=1:num_mc
    tmp_da_resp_mag = conv( hexa_data_an.da_resp_all.r(randperm(numel(hexa_data_an.da_resp_all.r)))  , trial_kernel , 'same' );
    tmp_cross_corr(kk,:) = xcorr(rew_rate-mean(rew_rate),tmp_da_resp_mag-mean(tmp_da_resp_mag),shifts );
end

shadedErrorBar(-shifts:shifts,mean(tmp_cross_corr,1),2.*std(tmp_cross_corr,[],1)); hold on;
plot(-shifts:shifts,cross_corr,'color',cat_map(1,:),'linewidth',4);
ylabel('Cross-correlation DA response and reward rate');
xlabel('Shift (trials)');

%% Some explicit dopamine analyses

cat_map = TNC_CreateRBColormap(8,'cat2');
cat_map = cat_map([1 3:6],:);

interportdist =      ...
[0	14	18	70	72.2	65.5 ;  ...
14	0	22.8	56	65.5	42; ...
18	22.8	0	72.2	70	56; ...
70	56	72.2	0	18	22.8;   ...
72.2	65.5	70	18	0	14; ...
65.5	42	56	22.8	14	0]

trans_p_dist.p = [];
trans_p_dist.d = [];
trans_p_dist.s = [];
Nback=1;

path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';
filenames = {'6PG5_NAc_conc_beh.csv'};


for session=1:2
    
    [hexa_data]     = HX_load_csv([path filenames{1}], 0, 1);
    [hexa_data_an]  = HX_analyze_session(hexa_data,session,1);
    
    tmp = find(sum(hexa_data_an.visits,1)==1);
    [~,visit_list_data] = max(hexa_data_an.visits(:,tmp(round(end/2):end)),[],1);
    [trans_mat_data] = HX_ComputeTransitionMatrix(visit_list_data,26,Nback);
    title(['DATA; Nback=' num2str(Nback) ' trans. matrix']);

    % create a basline probability for visitng each port
    for qq=1:6
        this_sess_base_p(qq) = numel(find(hexa_data_an.all_vis_ports==qq));
    end

    this_sess_base_p = this_sess_base_p./sum(this_sess_base_p); % normalize to P

    trans_p_dist.base_p(session,:) = this_sess_base_p;

    for qq=1:6
        for pp=1:6
            trans_mat_base(qq,pp) = this_sess_base_p(qq)*this_sess_base_p(pp);
        end
    end
    
    dist_lin = reshape(interportdist,1,36);
    trans_lin = reshape((trans_mat_data./sum(sum(trans_mat_data)))-trans_mat_base,1,36);
    
    trans_p_dist.p = [trans_p_dist.p trans_lin(dist_lin>0)];
    trans_p_dist.d = [trans_p_dist.d dist_lin(dist_lin>0)];
    trans_p_dist.s = [trans_p_dist.s ones(1,30)*session];

end

trans_p_dist.d_u = unique(trans_p_dist.d);

for mm=trans_p_dist.d_u
   trans_p_dist.p_avg(find(trans_p_dist.d_u==mm)) = mean(trans_p_dist.p(trans_p_dist.d==mm));
   trans_p_dist.p_err(find(trans_p_dist.d_u==mm)) = std(trans_p_dist.p(trans_p_dist.d==mm)) ./ sqrt(numel(find(trans_p_dist.d==mm)));
end
figure(100); clf;
% scatter(trans_p_dist.d,trans_p_dist.p,25,'k','filled','MarkerFaceAlpha',0.5);
% hold on;
errorbar(trans_p_dist.d_u,trans_p_dist.p_avg,trans_p_dist.p_err,'k','linewidth',2);
ylabel('P(visit, distance) - P(visit, previous)'); xlabel('Distance from previous port');

figure(101); clf;
for qq=1:5
    plot(1:6,trans_p_dist.base_p(qq,:),'color',cat_map(qq,:),'LineWidth',2); hold on;
end
box off; axis([0.9 6.1 0 0.6]);

% 1. if dopamine reflects knowledge of port quality then the 'inhibition' on un-rewarded visits should evolve to reflect rank ordering
    % would be surprising given that inhbition tends to develop more slowly
    % than positive RPE in Pavlovian paradigms
% 2. if there is a confidence like component - positive reward response magnitude should be inversely correlated to time since last visit/reward
% alternatively, could be more like a surprise / uncertainty effect so short intervals could yield increased magnitude

%
%--- so, need magnitude of da response in same space as visits(t,port)
%
