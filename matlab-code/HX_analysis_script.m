
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

clear model_compare

all_files = dir('~/Dropbox (HHMI)/hexaport/photometry/full_dataset/*conc_b*');
path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';

cost_per_port =                 ...
[0	14	18	70	72.2	65.5;   ...
14	0	22.8	56	65.5	42; ...
18	22.8	0	72.2	70	56; ...
70	56	72.2	0	18	22.8;   ...
72.2	65.5	70	18	0	14; ...
65.5	42	56	22.8	14	0]+0.1;

% cost_per_port = ones(6,6);

belief_model = 'p_check_match_win';
policy_model = 'e-proportional';
notes = 'da_compare_sandbox';
dir_path = [notes '_' belief_model '_' policy_model '/']
[SUCCESS,~,~] = mkdir(path,dir_path);

photo_flag = 1;

testing = 1;
for mmm = 3 %1:numel(all_files)
    
    breaks = strfind(all_files(mmm).name,'_');
    mouse_name = all_files(mmm).name(1:breaks(1)-1)
    
    session         = [1]; % session = 1;    
    photo_filename  = [path mouse_name '_photo.csv'];
    [hexa_data]     = HX_load_csv([path all_files(mmm).name], 0, photo_flag, photo_filename);
    [hexa_data_an]  = HX_analyze_session(hexa_data,session,photo_flag);

    % interim visualization
    figure(499); clf;
    subplot(121);
    port_color_map      = TNC_CreateRBColormap(8,'mapb');
    plot(hexa_data.event_time_con,[0 diff(hexa_data.session_n)'],'color',[0 0 0 0.25],'linewidth',2); hold on;
    for pp=1:6
        plot(hexa_data.event_time_con(hexa_data_an.visit_indices),hexa_data_an.p_choice_all(pp,:),'linewidth',2,'color',port_color_map(pp,:)); hold on;
    end
    axis([0 max(hexa_data.event_time_con(hexa_data_an.visit_indices)) 0 1]);
    box off;
    ylabel('P(visit,port)'); 

    model_compare.anim(mmm).mouse_name = mouse_name;
    model_compare.anim(mmm).belief_model = belief_model;
    model_compare.anim(mmm).policy_model = policy_model;
        
    Nback               = 1;

    intervals = [30 60 240 1200 2400];
    port_intervals = zeros(numel(session),6);
    for ss=session
        for qq=1:6
            port_intervals(ss,qq) = intervals(unique(hexa_data.port_rank(hexa_data.port_n==qq & ismember(hexa_data.session_n,ss))));
        end
        tmp                 = find(sum(hexa_data_an.visits,1)==1 & hexa_data_an.sessID'==ss);
        [~,visit_list_data] = max(hexa_data_an.visits(:,tmp),[],1);    
        subplot(1,numel(session),ss);
        [trans_mat_data]    = HX_ComputeTransitionMatrix(visit_list_data,0,Nback);
        exag = TNC_CreateRBColormap(8,'exag');
        imagesc(trans_mat_data,[0 0.2]); colormap(exag); axis equal; box off; colorbar;
        title(['DATA; Nback=' num2str(Nback) '; session: ' num2str(ss)]);
        [~,best_port]               = find(port_intervals==30);
        p_best_best                 = trans_mat_data(best_port(ss),best_port(ss));
        p_rew_best                  = sum(hexa_data_an.rewards(best_port(ss),:))./sum(hexa_data_an.visits(best_port(ss),:));

        model_compare.anim(mmm).pbb(ss) = p_best_best;
        model_compare.anim(mmm).base_stay(ss) = p_best_best/p_rew_best;
        model_compare.anim(mmm).trans(ss).trans_mat_data = trans_mat_data;
    end

    hist_wins = 500; %[0.1 1 2 10 20 50 100 500];

    for hh=1:numel(hist_wins)
        for reps=1:20
    
            [hexa_model]    = HX_model_session_2(hexa_data_an,policy_model,belief_model,cost_per_port,port_intervals,1,mean(model_compare.anim(mmm).base_stay),60*hist_wins(hh),0);  
    
            [port,event_time]   = find(hexa_model.visits==1);
            rewarded            = sum(hexa_model.rewards(:,event_time),1)';
            T                   = table(port,rewarded,event_time);
            writetable(T, [path dir_path 'Model_' mouse_name '_' policy_model '_' belief_model '_' num2str(reps) '.csv']);
    
            for ss=session
                % compare model predictions to data with transition matrix        
                figure(500);
                subplot(1,numel(session),ss);
        tmp                 = find(sum(hexa_data_an.visits,1)==1 & hexa_data_an.sessID'==ss);
        [~,visit_list_model] = max(hexa_model.visits(:,tmp),[],1);    

                % tmp = find(sum(hexa_data_an.visits,1)==1);
                % [~,visit_list_model] = max(hexa_model.visits(:,tmp),[],1);
                [trans_mat_model] = HX_ComputeTransitionMatrix(visit_list_model,0,Nback);
                imagesc(trans_mat_model,[0 0.2]); colormap(exag); axis equal; box off; colorbar;
                title(['MODEL; Nback=' num2str(Nback) ' Rho: ' num2str(corr2(trans_mat_data,trans_mat_model))]);
            end

            if reps==1
                tmp_mean_prew                               = zeros(size(hexa_model.p_reward,1),size(hexa_model.p_reward,2),10);
            end
            tmp_mean_trans(:,:,reps)                        = trans_mat_model;
            tmp_mean_prew(:,:,reps)                         = hexa_model.p_reward;
    
            model_compare.anim(mmm).corr(reps)              = corr2(trans_mat_data,trans_mat_model);
    
            model_compare.anim(mmm).rew_rate_model(reps)    = 100*(sum(sum(hexa_model.rewards)))/(sum(sum(hexa_model.rew_sched)));
            model_compare.anim(mmm).rew_rate_ideal(reps)    = 100*(sum(sum(hexa_model.ideal)))/(sum(sum(hexa_model.rew_sched)));
            model_compare.anim(mmm).rew_rate_random(reps)   = 100*(sum(sum(hexa_model.random)))/(sum(sum(hexa_model.rew_sched)));
        end
    
        model_compare.anim(mmm).rew_rate_mouse          = 100*(sum(sum(hexa_data_an.rewards)))/(sum(sum(hexa_model.rew_sched)));  
        
        % compare mean p_reward from model to DA responses observed in data
        all_rew     = find(hexa_data_an.da_visit_rew==1);
        all_rew_ids = hexa_data_an.da_visit_ids(all_rew);
        
        figure(4); clf;
        subplot(131);
            exag = TNC_CreateRBColormap(8,'exag');
            imagesc(mean(tmp_mean_trans,3),[0 0.2]); colormap(exag); axis equal; box off; colorbar;
        subplot(132);
            imagesc(mean(tmp_mean_prew,3));
        subplot(133);
            mn_p_rew = mean(tmp_mean_prew,3);
            visit_times = hexa_data.event_time(hexa_data_an.visit_indices);
            rew_times = visit_times(all_rew);
            % compute KL divergence between each neighboring point
            kl_div = zeros(1,size(mn_p_rew,2));
            for zz=2:size(mn_p_rew,2)
                kl_div(zz) = -sum( mn_p_rew(:,zz) .* log(mn_p_rew(:,zz)./mn_p_rew(:,zz-1)) );
            end
            % manhattan (cityblock) distance between neighboring p_reward estimates
            delta_p_rew = sum( abs( diff(mn_p_rew,2)),1); 
    
            [binnedData] = TNC_BinAndMean(delta_p_rew(hexa_data_an.visit_indices(all_rew)).*100, hexa_data_an.da_resp_all.r, 9);
            scatter(delta_p_rew(hexa_data_an.visit_indices(all_rew)),hexa_data_an.da_resp_all.r,25,[0.7 0.7 0.7],'filled'); hold on;
            errorbar(binnedData.bins.center./100,binnedData.bins.avg,binnedData.bins.sem,'color',[0 0 0 0.5],'LineWidth',3);
            axis([0 2 0 5]); xlabel('\Sigma |\DeltaP(rew|port)|'); ylabel('DA reward response');
        
            [rho,pval_rho] = corrcoef(delta_p_rew(hexa_data_an.visit_indices(all_rew)),hexa_data_an.da_resp_all.r);
            disp(' ');
            disp('_______________________________________');
            disp(['History length: ' num2str(hist_wins(hh))]);
            disp(['DA corr: ' num2str(rho(1,2)) ' ; p=' num2str(pval_rho(1,2))]);        
            disp(['BEHAV r2: ' num2str( mean(model_compare.anim(mmm).corr.^2) )]);        


figure(499);
subplot(122);

plot(hexa_data.event_time_con,[0 diff(hexa_data.session_n)'],'color',[0 0 0 0.25],'linewidth',2); hold on;
for pp=1:6
    plot(1:size(hexa_model.visits,2),cumsum(hexa_model.visits(pp,:))./cumsum(sum(hexa_model.visits,1)),'linewidth',2,'color',port_color_map(pp,:)); hold on;
end
axis([0 max(hexa_data.event_time_con(hexa_data_an.visit_indices)) 0 1]);
box off;
ylabel('P(visit,port)');             
title(['History length: ' num2str(hist_wins(hh))]);   

    end

    model_compare.anim(mmm).mean_trans = mean(tmp_mean_trans,3);
    model_compare.anim(mmm).mean_p_rew = mean(tmp_mean_prew,3);
    model_compare.anim(mmm).deltaprew  = delta_p_rew(hexa_data_an.visit_indices(all_rew));
    model_compare.anim(mmm).da_resp_m  = hexa_data_an.da_resp_all.r;
    model_compare.anim(mmm).da_resp_i  = all_rew_ids;    

end

if testing==0
    save([path dir_path 'model_compare_out'],'model_compare')
end

%% Looking at properties of different bandwidths of DA activity

    photo_sess_data     = sgolayfilt(hexa_data.photo.dFF(ismember(hexa_data.photo.sess,session)),3,51);
    port_color_map      = TNC_CreateRBColormap(8,'mapb');

    sym_map             = TNC_CreateRBColormap(1000,'bb-sym');
    
    photo_event_win     = [250 500];

    visit_indices       = find(hexa_data.unique_vis==1 & ~isnan(hexa_data.photo_i_con) & ismember(hexa_data.session_n,session));
 
    % need a correction for the fact that photo_i resets each session
    % whereas other indexing operations do not
    
    photo_visit_inds    = hexa_data.photo_i_con(visit_indices);
    port_visit_ids      = hexa_data.port_n(visit_indices);
    rew_visit_ids       = hexa_data.rewarded(visit_indices);
    rew_rate_in_visits  = conv(rew_visit_ids,[0 ones(1,50) 0]/50,'same');

    % 
    % filter_yes=1;
    % 
    % if filter_yes
    %     [lowBandData,hiBandData] = TNC_FilterPhotoData(photo_sess_data(1:1e5),100,0,1,[1 1]);
    % end

    range = 1:numel(photo_sess_data);
    d1 = designfilt("lowpassiir",FilterOrder=12, ...
        HalfPowerFrequency=0.0001,DesignMethod="butter");
    y = filtfilt(d1,photo_sess_data(range));
    d2 = designfilt("highpassiir",FilterOrder=12, ...
        HalfPowerFrequency=0.0001,DesignMethod="butter");
    yH = filtfilt(d2,photo_sess_data(range));

    delta_visit = zeros(1,numel(y));
    delta_reward = zeros(1,numel(y));
    delta_visit(photo_visit_inds) = 1;
    delta_reward(photo_visit_inds(rew_visit_ids==1)) = 1;

    kernel = [0 ones(1,5e3) 0]./5;
    visit_rate = conv(delta_visit,kernel,'same');
    rew_rate = conv(delta_reward,kernel,'same');

    figure(1); clf;
    % plot(photo_sess_data(range));
    plot(yH,'LineWidth',1);
    hold on;
    plot(y,'LineWidth',2);
    plot(visit_rate-4);
    plot(rew_rate-4);
    plot(photo_visit_inds,-3*ones(1,numel(photo_visit_inds)),'k*');
    plot(photo_visit_inds(rew_visit_ids==1),-3*ones(1,numel(photo_visit_inds(rew_visit_ids==1))),'c*');

    maxlag=2e5;
    low_by_rew = xcorr(y,rew_rate,maxlag);
    low_by_vis = xcorr(y,visit_rate,maxlag);
    figure(2); clf;
    plot(-maxlag:maxlag,low_by_rew);
    hold on;
    plot(-maxlag:maxlag,low_by_vis);
    

    photo_event_win     = [2500 2500];
    hp_sink             = TNC_ExtTrigWins(yH,photo_visit_inds,photo_event_win);
    lp_sink             = TNC_ExtTrigWins(y,photo_visit_inds,photo_event_win);
    bb_sink             = TNC_ExtTrigWins(photo_sess_data,photo_visit_inds,photo_event_win);

    sym_map_r = TNC_CreateRBColormap(1024,'wred');
    sym_map_b = TNC_CreateRBColormap(1024,'wblue');
    figure(3); clf;
    imagesc(-hp_sink.wins([find(rew_visit_ids==1) ; find(rew_visit_ids==0)],:),[-5 5]); colormap([sym_map_r;sym_map_b]);

    figure(4); clf;
    for kk=1:6
        plot(hp_sink.range,mean(hp_sink.wins(find(rew_visit_ids==1 & port_visit_ids==kk),:),1),'LineWidth',1,'color',[0.5 0.5 0.5]);
        hold on;
        plot(hp_sink.range,mean(hp_sink.wins(find(rew_visit_ids==0 & port_visit_ids==kk),:),1),'color',[0.9 0.3 0.3]);
    end

    plot(hp_sink.range,mean(hp_sink.wins(find(rew_visit_ids==1),:),1),'k','LineWidth',3);
    hold on;
    plot(hp_sink.range,mean(hp_sink.wins(find(rew_visit_ids==0),:),1),'r','LineWidth',3);
    % plot(bb_sink.range,mean(bb_sink.wins(find(rew_visit_ids==1),:),1),'LineWidth',2);
    % plot(bb_sink.range,mean(bb_sink.wins(find(rew_visit_ids==0),:),1),'LineWidth',2);
    yyaxis right;
    plot(lp_sink.range,mean(lp_sink.wins(find(rew_visit_ids==1),:),1),'k');
    plot(lp_sink.range,mean(lp_sink.wins(find(rew_visit_ids==0),:),1),'r');
    
%% Plot the compare data
clear compiled_data prediction pbb

for qq=1:numel(model_compare.anim)

    compiled_data(qq,1) = model_compare.anim(qq).rew_rate_mouse;
    compiled_data(qq,2) = mean(model_compare.anim(qq).rew_rate_model);
    compiled_data(qq,3) = mean(model_compare.anim(qq).rew_rate_ideal);
    compiled_data(qq,4) = mean(model_compare.anim(qq).rew_rate_random);

    prediction(qq) = mean(model_compare.anim(qq).corr .^ 2);
    pbb(qq) = model_compare.anim(qq).pbb(end);

end

figure(500); clf;
subplot(121);
cmap = TNC_CreateRBColormap(size(compiled_data,1),'cpb');
% boxplot(( (compiled_data-compiled_data(:,4))./ (compiled_data(:,3)-compiled_data(:,4)) ),{'Data' 'Model' 'Ideal' 'Random'}); axis([0 5 -0.2 1]); hold on;
boxplot(compiled_data,{'Data' 'Model' 'Ideal' 'Random'}); axis([0 5 0 100]); hold on;
for zz=1:size(compiled_data,1)
    plot(1:4,compiled_data(zz,1:4),'color',cmap(zz,:))
end
ylabel('Collection efficiency (optimality index)'); box off;

subplot(122);
boxplot(prediction',{'Transition matrix'}); axis([0 3 0 1]);
ylabel('R^2'); box off;

figure(501); clf;
scatter(prediction,pbb,25,'k','filled');
xlabel('R2 model prediction'); ylabel('P(stay at best port)');
axis([0 1 0 0.5])

[P,ANOVATAB,STATS] = anova1(compiled_data,{'Data' 'Model' 'Ideal' 'Random'});
COMPARISON = multcompare(STATS,'alpha',0.01,"CriticalValueType","dunnett")

%% Plotting output of most recent model run

run_new = 1;
all_files = dir('~/Dropbox (HHMI)/hexaport/photometry/full_dataset/*conc_b*');
path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';

if run_new
    mmm = 2 %randperm(numel(all_files),1)

    cost_per_port =      ...
        [0	14	18	70	72.2	65.5 ;  ...
        14	0	22.8	56	65.5	42; ...
        18	22.8	0	72.2	70	56; ...
        70	56	72.2	0	18	22.8;   ...
        72.2	65.5	70	18	0	14; ...
        65.5	42	56	22.8	14	0].^1.5;

    session = [1 2]; % session = 1;   
    belief_model = 'p_check_match';
    policy_model = 'e-proportional';

    % belief_model = 'p_check_match';
    % policy_model = 'softmax';
    
    [hexa_data]     = HX_load_csv([path all_files(mmm).name], 0, 1);
    % figure out port ranks for model and analysis
    intervals = [30 60 240 1200 2400];
    port_intervals = zeros(1,6);
    for qq=1:6
        port_intervals(qq) = intervals(unique(hexa_data.port_rank(hexa_data.port_n==qq & ismember(hexa_data.session_n,session))));
    end

    [hexa_data_an]  = HX_analyze_session(hexa_data,session,0);
    [hexa_model]    = HX_model_session_2(hexa_data_an,policy_model,belief_model,cost_per_port,port_intervals,1,0);  
end

%/---------- DATA plotting
figure(100); clf;
cat_map = TNC_CreateRBColormap(8,'mapb');
all_vis_rate = [0 diff(sgolayfilt(cumsum(sum(hexa_data_an.visits,1)),3,1001))];

subplot(122);
plot([0 1],[0 1],'k'); hold on;

visit_win   = 2500;
all_visits  = conv(sum(hexa_data_an.visits,1),[0 ones(1,visit_win) 0],'same');
all_rews    = conv(sum(hexa_data_an.rewards,1),[0 ones(1,visit_win) 0],'same');

for qq=1:6 


    this_port_visits        = conv(hexa_data_an.visits(qq,:),[0 ones(1,visit_win) 0],'same');
    this_port_rews          = conv(hexa_data_an.rewards(qq,:),[0 ones(1,visit_win) 0],'same');

    frac_visits_this_port   = this_port_visits./all_visits;
    frac_rew_this_port      = this_port_rews./all_rews;

    % frac_visits_this_port   = cumsum(hexa_data_an.visits(qq,:))./cumsum(sum(hexa_data_an.visits,1));
    % frac_rew_this_port      = cumsum(hexa_data_an.rewards(qq,:))./cumsum(sum(hexa_data_an.rewards,1));

    subplot(121);
    plot(1:size(hexa_data_an.visits(qq,:),2),frac_visits_this_port,'color',[cat_map(qq,:)],'LineWidth',3); hold on; box off;

    subplot(122);
    scatter(frac_rew_this_port,frac_visits_this_port,50,numel(frac_rew_this_port):-1:1,'filled','MarkerEdgeColor',cat_map(qq,:)); colormap(bone); hold on;
    axis([0 0.7 0 0.7]); xlabel('Fraction rewards'); ylabel('Fraction visits');
end

subplot(121);
plot(1:size(hexa_data_an.visits(qq,:),2),[0 diff(hexa_data_an.sessID')],'k');
axis([0 numel(frac_rew_this_port) 0 0.7]); ylabel('Fraction visits'); xlabel('Time');

%/---------- MODEL plotting
figure(101); clf;
cat_map = TNC_CreateRBColormap(8,'mapb');
all_vis_rate = [0 diff(sgolayfilt(cumsum(sum(hexa_model.visits,1)),3,1001))];

subplot(122);
plot([0 1],[0 1],'k'); hold on;

visit_win   = 2500;
all_visits  = conv(sum(hexa_model.visits,1),[0 ones(1,visit_win) 0],'same');
all_rews    = conv(sum(hexa_model.rewards,1),[0 ones(1,visit_win) 0],'same');

for qq=1:6 


    this_port_visits        = conv(hexa_model.visits(qq,:),[0 ones(1,visit_win) 0],'same');
    this_port_rews          = conv(hexa_model.rewards(qq,:),[0 ones(1,visit_win) 0],'same');

    frac_visits_this_port   = this_port_visits./all_visits;
    frac_rew_this_port      = this_port_rews./all_rews;

    % frac_visits_this_port   = cumsum(hexa_data_an.visits(qq,:))./cumsum(sum(hexa_data_an.visits,1));
    % frac_rew_this_port      = cumsum(hexa_data_an.rewards(qq,:))./cumsum(sum(hexa_data_an.rewards,1));

    subplot(121);
    plot(1:size(hexa_model.visits(qq,:),2),frac_visits_this_port,'color',[cat_map(qq,:)],'LineWidth',3); hold on; box off;

    subplot(122);
    scatter(frac_rew_this_port,frac_visits_this_port,50,numel(frac_rew_this_port):-1:1,'filled','MarkerEdgeColor',cat_map(qq,:)); colormap(bone); hold on;
    axis([0 0.7 0 0.7]); xlabel('Fraction rewards'); ylabel('Fraction visits');
end

subplot(121);
plot(1:size(hexa_model.visits(qq,:),2),[0 diff(hexa_data_an.sessID')],'k');
axis([0 numel(frac_rew_this_port) 0 0.7]); ylabel('Fraction visits'); xlabel('Time');


cat_map = TNC_CreateRBColormap(8,'cat2');
cat_map = cat_map([1 3:6 8],:);
trial_win = 2;
trial_kernel = TNC_CreateGaussian(500,trial_win,1000,1);

Nback=1;

tmp = find(sum(hexa_data_an.visits,1)==1);
[visit_times,visit_list_data] = max(hexa_data_an.visits(:,tmp),[],1);
[trans_mat_data] = HX_ComputeTransitionMatrix(visit_list_data,26,Nback);
title(['DATA; Nback=' num2str(Nback) ' trans. matrix']);

tmp = find(sum(hexa_data_an.visits,1)==1);
[~,visit_list_model] = max(hexa_model.visits(:,tmp),[],1);
[trans_mat_model] = HX_ComputeTransitionMatrix(visit_list_model,27,Nback);
title(['MODEL; Nback=' num2str(Nback) ' Rho: ' num2str(corr2(trans_mat_data,trans_mat_model))]);


exag = TNC_CreateRBColormap(8,'exag');
figure(28); imagesc(cost_per_port,[0 75]); colormap(exag);

figure(30); clf;
% plot the evolution of P(1_n,1_n-1) transition
visit_list = visit_list_data(1:end-1);
for zz=1
    this_port_visits = find(visit_list==zz);
    for mm=1
        p11_trans = this_port_visits( find(visit_list_data(this_port_visits+1)==mm) );
    end
end
pN1_visits = zeros(1,numel(visit_list_data));
pN1_visits(this_port_visits) = 1;
p11_visits = zeros(1,numel(visit_list_data));
p11_visits(p11_trans) = 1;

plot(conv(p11_visits,[0 ones(1,300) 0]/300,'same')./conv(pN1_visits,[0 ones(1,300) 0]/300,'same')); hold on;
title(['DATA; p(1,1) / p(N,1)']);  axis([0 numel(visit_list_model) 0 0.75]);

figure(31); clf;
% plot the evolution of P(1_n,1_n-1) transition
visit_list = visit_list_model(1:end-1);
for zz=1
    this_port_visits = find(visit_list==zz);
    for mm=1
        p11_trans = this_port_visits( find(visit_list_model(this_port_visits+1)==mm) );
    end
end
pN1_visits = zeros(1,numel(visit_list_model));
pN1_visits(this_port_visits) = 1;
p11_visits = zeros(1,numel(visit_list_model));
p11_visits(p11_trans) = 1;

plot(conv(p11_visits,[0 ones(1,300) 0]/300,'same')./conv(pN1_visits,[0 ones(1,300) 0]/300,'same')); hold on;
title(['MODEL; p(1,1) / p(N,1)']); axis([0 numel(visit_list_model) 0 0.75]);

figure(40); clf;
plot(2:6,trans_mat_data(2:6,1)); hold on;
plot(2:6,trans_mat_model(2:6,1)); hold on;
axis([2 6 0 0.25]); legend({'data','model'}); ylabel('P(N->1)'); xlabel('N');

% compare the probability of reward on visits to best port coming from same
% port vs any port
% for port_id = 1:3
%     port_id
%     stay_rewards = sum(hexa_model.rewards(port_id,hexa_model.visits(port_id,:)==1 & hexa_model.stay_go==1)) ./ numel(find(hexa_model.visits(port_id,:)==1 & hexa_model.stay_go==1))
%     return_rewards = sum(hexa_model.rewards(port_id,hexa_model.visits(port_id,:)==1 & hexa_model.stay_go==0)) ./ numel(find(hexa_model.visits(port_id,:)==1 & hexa_model.stay_go==0))
%     hexa_model.p_reward(port_id,end-1)
% end

%% Plotting output including photometry from above

figure(200); clf;
visit_indices       = find(hexa_data.unique_vis==1 & ~isnan(hexa_data.photo_i) & ismember(hexa_data.session_n,session));
rew_visit_ids       = hexa_data.rewarded(visit_indices);

scatter(hexa_data_an.da_resp_all.t,hexa_data_an.da_resp_all.r,50,hexa_data_an.da_resp_all.p,'filled','MarkerFaceAlpha',0.25); colormap(cat_map); hold on;
plot(hexa_data_an.da_resp_all.t,conv( hexa_data_an.da_resp_all.r  , trial_kernel , 'same' ) ,'color' , [0 0.67 1 0.5], 'linewidth' , 4);
axis([0 max(visit_indices) -500 1500]);

yyaxis right;
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

%% Code to convert model output into Laura's table version

[port,event_time] = find(hexa_model.visits==1);
rewarded = sum(hexa_model.rewards(:,event_time),1)';
T=table(port,rewarded,event_time);
writetable(T, 'ModelTestExport_6PG12.csv');

