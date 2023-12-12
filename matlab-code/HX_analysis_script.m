
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
% [~,visit_list] = max(hexa_model.visits(:,tmp),[],1);
[~,visit_list] = max(hexa_model.visits(:,tmp(1:round(end/3))),[],1);
% [~,visit_list] = max(hexa_model.visits(:,tmp(round(end/3):end)),[],1);


[trans_mat] = HX_ComputeTransitionMatrix(visit_list,25,Nback);
title(['SIMULATION; Nback=' num2str(Nback) ' trans. matrix']);


tmp = find(sum(hexa_data_an.visits,1)==1);
% [~,visit_list_data] = max(hexa_data_an.visits(:,tmp),[],1);
[~,visit_list_data] = max(hexa_data_an.visits(:,tmp(1:round(end/3))),[],1);
% [~,visit_list_data] = max(hexa_data_an.visits(:,tmp(round(end/3):end)),[],1);

[trans_mat_data] = HX_ComputeTransitionMatrix(visit_list_data,26,Nback);
title(['DATA; Nback=' num2str(Nback) ' trans. matrix']);

