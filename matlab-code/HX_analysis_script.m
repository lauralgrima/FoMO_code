
%% Requirements for figure
% Model schematic: 
    % Possible policies p(sample) and p(trans) - choice based on these (illustrate the differences between them)

% 
% Paired comparison of reward rate for p(sample) and p(trans) after random - one dot per mouse 
% 
% Delta reward rate (mouse - model) for p(sample and p(trans) - might be combined with b?
% 
% Reward x visit ratio matching plot (maybe one each for p(sample) and p(trans)?) 
% and then sensitivity calculation to show undermatching - maybe compared to empirical sensitivity which should be very similar 
% 
% Sensitivity slope values split for conditional matching, model vs. animal
% 
% R2 for transition matrix 

%% Summary statistics for best fit model

% choice ratio is choices to port N divided by choices to port ~N (log10)

% Load each model run, find best option 

% calculate matching statistics and cumulative choice plots to compare to
% data

% optimal R2

% income RMSE

%% Create analyzed versions of all data files for subsequent use and fitting

clear model_compare

all_files = dir('~/Dropbox (HHMI)/hexaport/photometry/full_dataset/*conc_b*');

% example special case if want to run just on single animal
% all_files = dir('~/Dropbox (HHMI)/hexaport/photometry/full_dataset/6PG31_NAc_conc_b*');

path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';

cost_per_port =                 ...
[1	14	18	70	72.2	65.5;   ...
14	1	22.8	56	65.5	42; ...
18	22.8	1	72.2	70	56; ...
70	56	72.2	1	18	22.8;   ...
72.2	65.5	70	18	1	14; ...
65.5	42	56	22.8	14	1];

session         = 1
notes = ['da_store_analyzed_sess' num2str(session) 'nLL_RewCnt_sess3'];
dir_path = [notes '/']
[SUCCESS,~,~] = mkdir(path,dir_path);

photo_flag = 1;
figure(600);

port_color_map      = TNC_CreateRBColormap(8,'mapb');

testing = 1;
for mmm = 1:numel(all_files) % mice 11 and 16 do not have session 2 data
    
    clear prior;
    breaks = strfind(all_files(mmm).name,'_');
    mouse_name = all_files(mmm).name(1:breaks(1)-1);
      
    photo_filename  = [path mouse_name '_photo.csv'];
    [hexa_data]     = HX_load_csv([path all_files(mmm).name], 0, photo_flag, photo_filename);

    if min(ismember(session,unique(hexa_data.session_n)))==1

            [hexa_data_an]  = HX_analyze_session(hexa_data,session,photo_flag);
        
            intervals = [30 60 240 1200 2400];
            port_intervals = zeros(numel(session),6);
            for ss=session
                for qq=1:6
                    port_intervals(ss,qq)       = intervals(unique(hexa_data.port_rank(hexa_data.port_n==qq & ismember(hexa_data.session_n,ss))));
                    port_rank_this_sess(ss,qq)  = (unique(hexa_data.port_rank(hexa_data.port_n==qq & ismember(hexa_data.session_n,ss))));
                end
            end    
        
        %----------------------------------------
        %----------------------------------------
        %-------------- FITTING ALPHA PER SESSION
        %----------------------------------------
        %----------------------------------------

        visit_matrix = hexa_data_an.visits;

        sample_logic = sum(visit_matrix,1);
        
        all_visits = find(sum(hexa_data_an.visits,1)==1);
        rew_logic = sum(hexa_data_an.rewards,1);
        all_rewards = rew_logic(all_visits);
        income = movmean(all_rewards,51);
        
        % sampling rate is now set to 1 Hz
        frame_rate = 1;
        
        hexa_model.rew_sched = zeros(size(hexa_data_an.visits));
        for ss=unique(hexa_data_an.sessID)'
            valid_inds = find(hexa_data_an.sessID==ss);
            for qq=1:6
                hexa_model.rew_sched(qq,valid_inds(1):round(port_intervals(ss,qq)*frame_rate):valid_inds(end)) = 1;
            end
        end
        hexa_model.rew_sched(:,2) = 1;
        
        rew_sched = hexa_model.rew_sched;
        
        % example values for alpha_params
        alpha_params_init = [0 0.2 0.2 300 500];
        v_ind = 1:sum(sample_logic);

        vis_inds                = find(sum(visit_matrix,1)==1);
        [~,visit_list_data]     = max(visit_matrix(:,vis_inds),[],1);    
        [trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data(1:end),0,1);
        
        % Visualize the alpha function with default
        % alpha_version = 'doub_exp'
        % alpha_vis_init = alpha_params_init(1) + (alpha_params_init(2)*(1-exp(-v_ind/alpha_params_init(4))) .* (alpha_params_init(3)*exp(-v_ind/alpha_params_init(5))));
        
        alpha_version = 'sig_exp';
        [rise_kern] = TNC_CreateGaussian(alpha_params_init(4),100,sum(sample_logic),1);
        alpha_rise = cumsum(rise_kern)*alpha_params_init(2);
        alpha_vis_init = alpha_params_init(1) + (alpha_rise .* (alpha_params_init(3)*exp(-v_ind/alpha_params_init(5))));
        
        % Online plot initialization
        exag_map    = TNC_CreateRBColormap(8,'exag');
        figure(10); clf; 
        figure(11); clf;
        
        % Examine the base r2 with a guess at alpha initialization
        a1 = alpha_params_init(1);
        a2 = alpha_params_init(2);
        a3 = alpha_params_init(3);
        a4 = alpha_params_init(4);
        a5 = alpha_params_init(5);
        
        num_iter        = 8;
        vismat          = zeros(6,numel(all_visits),num_iter);
        rewmat          = zeros(6,numel(all_visits),num_iter);
        trans_mat_modrun= zeros(6,6,num_iter);

        if session==1
            % assume random initialization
            prior(:,1) = 1/6;
            prior(:,2) = 0.1;
        else
            % estimate prior from hexa_data_an for session-1
            [hexa_data_an_prior]    = HX_analyze_session(hexa_data,session-1,0);
            tmp                     = find(sum(hexa_data_an_prior.visits,1)==1);
            [~,visit_list_data]     = max(hexa_data_an_prior.visits(:,tmp(round(0.7*end:end))),[],1);        
            [trans_mat_data_prior]  = HX_ComputeTransitionMatrix(visit_list_data,0,1);
            prior(:,2)              = diag(trans_mat_data_prior,0);
            prior(:,1)              = sum(trans_mat_data_prior,1)' - diag(trans_mat_data_prior,0);
            % add noise?
        end


        %----------------------------------------
        %-------------- FITTING ALPHA TO DA
        %----------------------------------------
        alpha = @(a1,a2,a3,a4,a5,x) a1 + (a2 ./ (1+exp((a4-x)/(a4./6)))) .*  (a3*exp(-x/a5));
        fitfun = fittype( alpha );        
        targety = movmean(hexa_data_an.da_resp_all.r,3)./max(movmean(hexa_data_an.da_resp_all.r,11));
        a0 = [ 0 0.5 0.5 100 1000 ];        
        [dopa.f,dopa.gof] = fit([1:numel(hexa_data_an.da_resp_all.r)]',targety,fitfun,'StartPoint',a0,'Upper',[0.1 1 1 numel(targety) 2*numel(targety)],'Lower',[0 0 0 20 20]);

        %----------------------------------------
        %-------------- FITTING ALPHA TO DA
        %----------------------------------------
        % grid search optimization
        a2_vec = [0.05 0.1 0.2 0.5 0.99];
        a4_vec = [10 20 50 100 200];
        a5_vec = [50 100 200 500 1000];

        [~,close_a4] = min(abs(dopa.f.a4-a4_vec));
        [~,close_a5] = min(abs(dopa.f.a5-a5_vec));

        for a2 = a2_vec
            for a4 = a4_vec
                for a5 = a5_vec
        
                    a1 = 0.001;
                    a2 = a2;
                    a3 = a2;
                    trans_r2_iter  = zeros(1,num_iter);
                    income_r2_iter = zeros(1,num_iter);

                    parfor iter = 1:num_iter
                        [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(a1,a2,a3,a4,a5,alpha_version,visit_matrix,cost_per_port,rew_sched,income,prior);
                    end
        
                    %--------
                    % Maybe just update the fit to -LL of observed choices using n_iters to get an estimate of P(obs choice | model)?
                    % I think that would just be log( observed visit matrix - estimated
                    % probability of choices from simulations ) -> summed over total visits or
                    % mean per visit
                    opt_r2_tensor(find(a2==a2_vec),find(a4==a4_vec),find(a5==a5_vec)) = median(trans_r2_iter);
                    opt_inc_tensor(find(a2==a2_vec),find(a4==a4_vec),find(a5==a5_vec)) = median(income_r2_iter);
                    params_a2(find(a2==a2_vec),find(a4==a4_vec),find(a5==a5_vec)) = a2;
                    params_a5(find(a2==a2_vec),find(a4==a4_vec),find(a5==a5_vec)) = a5;

                    % compare observed visits to model expectations
                    % Switch to using KL divergence (which is LL estimator) 
                    % old nLL version
                    % vis_prob    = mean(vismat,3);                    
                    % vis_obs     = hexa_data_an.visits(:,all_visits);
                    % nLL         = -mean( log( vis_prob(vis_obs==1)+0.001 ) );

                    LL = sum( (reshape(trans_mat_data+0.001,1,36)) .* log( reshape(trans_mat_data+0.001,1,36) ./ reshape(squeeze(mean(trans_mat_modrun,3)+0.001),1,36) ) );
                    opt_LL_tensor(find(a2==a2_vec),find(a4==a4_vec),find(a5==a5_vec)) = mean(LL);

                    tot_rew     = sum(sum(mean(rewmat,3)));               
                    opt_RColl_tensor(find(a2==a2_vec),find(a4==a4_vec),find(a5==a5_vec)) = tot_rew;
                    
                    % alpha_vis = a1 + (a2*(1-exp(-v_ind/a4)) .* (a3*exp(-v_ind/a5)));
                    % figure(10); subplot(1,numel(a2_vec),find(a2==a2_vec)); plot(v_ind,alpha_vis); hold on; axis([v_ind(1) v_ind(end) 0 max(a2_vec).^2]); box off;
                    if find(a5==a5_vec)==numel(a5_vec)
                        hhh = figure(11); 
                        subplot(4,numel(a2_vec),find(a2==a2_vec));
                        imagesc(squeeze(opt_r2_tensor(find(a2==a2_vec),:,:)),[0.25 0.85]); colormap(exag_map);                        
                        title('Trans R2');
                        figure(11); hold off; subplot(4,numel(a2_vec),find(a2==a2_vec)+numel(a2_vec));
                        imagesc(squeeze(opt_inc_tensor(find(a2==a2_vec),:,:)),[0.05 0.25]); colormap(exag_map);
                        hold on;
                        scatter(close_a5,close_a4,50,'k','filled');
                        title('Income RMSE');
                        figure(11); subplot(4,numel(a2_vec),find(a2==a2_vec)+numel(a2_vec)+numel(a2_vec));
                        imagesc(squeeze(opt_LL_tensor(find(a2==a2_vec),:,:)),[0.5 5]); colormap(exag_map);
                        title('nLL');
                        figure(11); subplot(4,numel(a2_vec),find(a2==a2_vec)+numel(a2_vec)+numel(a2_vec)+numel(a2_vec));
                        imagesc(abs(squeeze(opt_RColl_tensor(find(a2==a2_vec),:,:))-sum(sum(hexa_data_an.rewards))),[0 500]); colormap(exag_map);
                        title(['Total Reward - mouse: ' num2str(sum(sum(hexa_data_an.rewards)))]);
                    end
                end
            end
        end

        exportgraphics(hhh, [path dir_path mouse_name '_FitTensSummary.pdf'],"ContentType","vector"); % write out Fig 11
        
        % look for joint min
        summary_fit_fig = figure(700);
        subplot(ceil(numel(all_files)/5),5,mmm)
        scatter(reshape(opt_inc_tensor,1,5*5*5),reshape(opt_r2_tensor,1,5*5*5),(1+reshape(params_a2,1,5*5*5)).^2*50,reshape(params_a5,1,5*5*5),'filled','MarkerEdgeColor','k'); colormap(exag_map);
        xlabel('inc RMSE'); ylabel('Transition matrix r^2');
        title([mouse_name '; r2: ' num2str(max(reshape(opt_r2_tensor,1,5*5*5))) '; rmse: ' num2str(min(reshape(opt_LL_tensor,1,5*5*5)))]);
        axis([0.05 0.25 0 1]);
        
        %----------------------------------------
        %----------------------------------------
        %-------------- FITTING ALPHA PER SESSION
        %----------------------------------------
        %----------------------------------------
        
           % interim visualization
            Nback               = 1;
            hh = figure(499); clf;
            subplot(1,3,1:2);
            % plot(hexa_data.e√vent_time_con,[0 diff(hexa_data.session_n)'],'color',[0 0 0 0.25],'linewidth',2); hold on;
            for pp=1:6
                plot(hexa_data.event_time_con(hexa_data_an.visit_indices),hexa_data_an.p_choice_all(pp,:),'linewidth',2,'color',port_color_map(unique(hexa_data.port_rank(hexa_data.port_n==pp & ismember(hexa_data.session_n,ss))),:)); hold on;
            end
            axis([0 max(hexa_data.event_time_con(hexa_data_an.visit_indices)) 0 0.67]);
            box off;
            ylabel('P(visit,port)'); 
            title(mouse_name);
            
            [~,ii_r]            = sort(port_rank_this_sess(1,:));
            tmp                 = find(sum(hexa_data_an.visits,1)==1);
            [~,visit_list_data] = max(hexa_data_an.visits(ii_r,tmp),[],1);    
        
            [trans_mat_data]    = HX_ComputeTransitionMatrix(visit_list_data,0,Nback);
            exag = TNC_CreateRBColormap(8,'exag');
        
            subplot(133);
            imagesc(trans_mat_data,[0 0.2]); 
            colormap(exag); axis equal; box off; colorbar;
            title([mouse_name '; Nback=' num2str(Nback)]);
            drawnow;
        
            % exportgraphics(hh, [path dir_path mouse_name '_summary.pdf'],"ContentType","vector");
            % exportgraphics(hexa_data_an.da_hand1, [path dir_path mouse_name '_summary.pdf'],"ContentType","vector",'Append',true);
            % exportgraphics(hexa_data_an.da_hand2, [path dir_path mouse_name '_summary.pdf'],"ContentType","vector",'Append',true);
        
            % Proper per animal summary:
            % DA dynamics per trial, responses per port
            % Optimal R2 on trans matrix and income
            % Transition to suprarandom choices
        
            summary_fig = figure(600);    
            subplot(ceil(numel(all_files)/5),5,mmm)
            imagesc(trans_mat_data,[0 0.25]); colormap(exag); axis equal; box off; colorbar;
            title([ mouse_name '; r# ' num2str(sum(sum(hexa_data_an.rewards,1))) '; v# ' num2str(sum(sum(hexa_data_an.visits,1))) '; r2: ' num2str(max(reshape(opt_r2_tensor,1,5*5*5)))]);
            drawnow;
        
            eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_sess' num2str(session) '_opt.mat a* opt* dopa tot_rew visit_matrix cost_per_port rew_sched income prior port_rank_this_sess']);
            eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_sess' num2str(session) '_an.mat hexa_data_an']);
            disp(['Completed fitting for ' all_files(mmm).name ' session(s): ' num2str(session)]);

    else

        disp(['Skipped ' mouse_name ' because data lacked one of the target sessions; ' 'target session array: ' num2str(session) ' sessions found: ' num2str(unique(hexa_data.session_n)')])

    end

end

exportgraphics(summary_fig, [path dir_path 'all_mouse_summary_trans_mat.pdf'],"ContentType","vector");
exportgraphics(summary_fit_fig, [path dir_path 'all_mouse_summary_fit_mat.pdf'],"ContentType","vector");

%% Take saved files and calculate summary stats for paper

a=1;
b=0;
c=0;
frac = 0.95;
clear opt

sess=2;

figure(799+sess); clf;    
figure(9); clf;
% figure(850); clf;
figure(11); clf;

all_coms = [];
all_isos = [];
all_taus = [];
all_recloc = [];
all_r2 = [];
all_nLL = [];

all_sess_files = dir(['*sess' num2str(sess) '*_opt.mat']);
opt.r2      = zeros(numel(all_sess_files),9); % alpha fit, com fit, da fit, da raw, q, oio, wsls, random
opt.inc     = zeros(numel(all_sess_files),9); % alpha fit, com fit, da fit, da raw, q, oio, wsls, random
opt.rank    = zeros(numel(all_sess_files),1);
opt.rew     = zeros(numel(all_sess_files),1);

all_match       = zeros(6,2,numel(all_sess_files));
all_match_sens  = zeros(numel(all_sess_files),1);

for zz=1:numel(all_sess_files)

    clear S Z
    S = load(all_sess_files(zz).name);

    breaks = strfind(all_sess_files(zz).name,'_');

    loss = a*S.opt_r2_tensor-b*S.opt_inc_tensor-c*S.opt_LL_tensor;
    
    [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));
    sym = TNC_CreateRBColormap(numel(top_xperc_inds),'cpb');
        
    [a2_inds,a4_inds,a5_inds] = ind2sub(size(S.opt_r2_tensor),top_xperc_inds);
    
    [com_in_param_space] = centerOfMass3D(S.a2_vec(a2_inds), S.a4_vec(a4_inds), S.a5_vec(a5_inds), loss(top_xperc_inds)');

    % figure(900);
    % subplot(5,4,zz);
    % scatter(S.a4_vec(a4_inds), S.a5_vec(a5_inds),100,loss(top_xperc_inds)','filled'); colormap("parula");
    % xlim([min(S.a4_vec) max(S.a4_vec)]);
    % ylim([min(S.a5_vec) max(S.a5_vec)]);

    figure(799+sess);
    subplot(121);
    % scatter3(a4_vec(a4_inds),a5_vec(a5_inds),a2_vec(a2_inds),50,loss(top_xperc_inds),'filled'); colormap(sym);
    hold on;
    scatter3(com_in_param_space(2),com_in_param_space(3),com_in_param_space(1),100,mean(S.opt_r2_tensor(top_xperc_inds)'),'filled'); grid on;
    axis([min(S.a4_vec) max(S.a4_vec) min(S.a5_vec) max(S.a5_vec) min(S.a2_vec) max(S.a2_vec) ]); colormap(sym); colorbar;
    xlabel('a4'); ylabel('a5'); zlabel('a2|a3'); view([-40 15]);

    subplot(122);
    x = 1:1000;
    alpha_vis = 0.001 + (com_in_param_space(1) ./ (1+exp((com_in_param_space(2)-x)/(com_in_param_space(2)./6)))) .*  (com_in_param_space(1)*exp(-x/com_in_param_space(3)));    
    where_r2 = 2*(max(S.opt_r2_tensor(top_xperc_inds)')-0.5);
    plot(1:1000,alpha_vis,'color',[where_r2,0.5.*(1-where_r2),(1-where_r2)]); hold on;
    axis([0 1000 0 0.33]);
    ylabel('alpha(rew)'); xlabel('Num rewards'); box off;
    title([all_sess_files(zz).name(1:breaks(1)-1) ' ' all_sess_files(zz).name(breaks(2)-2:breaks(2)-1)]);
    drawnow; pause(0.1);

    all_r2      = [all_r2 ; max(S.opt_r2_tensor,[],"all")];
    all_coms    = [all_coms ; com_in_param_space];
    all_isos    = [all_isos ; numel(top_xperc_inds)];

    % load the corresponding hexa_data_an struct
        Z = load([all_sess_files(zz).name(1:end-7) 'an.mat']);

    % Get data/optimal discovered fit observations for comparison
    opt.r2(zz,1) = max(S.opt_r2_tensor,[],"all");
    opt.labels{1} = 'AQUA opt grid';
    opt.rew(zz,1) = S.tot_rew;
    opt.rew_act(zz,1) = numel(find(Z.hexa_data_an.rewards==1));


% ------------- ALPHA fitting routine

        alpha = @(a1,a2,a3,a4,a5,x) a1 + (a2 ./ (1+exp((a4-x)/(a4./6)))) .*  (a3*exp(-x/a5));
        fitfun = fittype( alpha );

        targety = movmean(Z.hexa_data_an.da_resp_all.r,3)./max(movmean(Z.hexa_data_an.da_resp_all.r,11));

        % reasonable initial guesses
        a0 = [ 0 0.5 0.5 100 1000 ];

        [f,gof] = fit([1:numel(Z.hexa_data_an.da_resp_all.r)]',targety,fitfun,'StartPoint',a0,'Upper',[0.1 1 1 numel(targety) 2*numel(targety)],'Lower',[0 0 0 10 20]);        

% ------------- ALPHA fitting routine

% ------------- 
% ------------- RUN RANGE OF MODEL SIMULATIONS

        num_iter = 24;

        all_visits = find(sum(Z.hexa_data_an.visits,1)==1);
        rew_logic = sum(Z.hexa_data_an.rewards,1);
        all_rewards = rew_logic(all_visits);

        vismat              = zeros(6,numel(all_visits),num_iter);
        rewmat              = zeros(6,numel(all_visits),num_iter);
        trans_r2_iter       = zeros(1,num_iter);
        income_r2_iter      = zeros(1,num_iter);

        [~,close_a4] = min(abs(f.a4-S.a4_vec))
        [~,close_a5] = min(abs(f.a5-S.a5_vec))

        prior = S.prior;
        
        [~,a2_ind_DA] = max(S.opt_r2_tensor(:,close_a4,close_a5));
        best_a2 = S.a2_vec(a2_ind_DA);

        % run AQUA for optimal alpha fit
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(0.001,com_in_param_space(1),com_in_param_space(1),com_in_param_space(2),com_in_param_space(3),'sig_exp',S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,2) = mean(trans_r2_iter);
        opt.rank(zz,1) = numel(find(S.opt_r2_tensor<=mean(trans_r2_iter)))./prod(size(S.opt_r2_tensor)); % percentile
        opt.labels{2} = 'AQUA opt com';
        opt.rew(zz,2) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,2) = mean(income_r2_iter);

        % range = [round(size(vismat,2)/2) : size(vismat,2)];
        range = [1 : size(vismat,2)];
        sum_vismat = sum(squeeze(mean(vismat(:,range,5),3)),2);
        sum_rewmat = sum(squeeze(mean(rewmat(:,range,5),3)),2);
        [~,port_ranki] = sort(sum(S.rew_sched,2),'descend');

        for jjj=1:6
            all_match(jjj,2,zz) = log10( sum_vismat(port_ranki(jjj)) ./ sum(sum_vismat(port_ranki~=port_ranki(jjj))) ); % choice_ratio 
            all_match(jjj,1,zz) = log10( sum_rewmat(port_ranki(jjj)) ./ sum(sum_rewmat(port_ranki~=port_ranki(jjj))) ); % reward_ratio 
        end

        [coeffs] = polyfit(all_match(:,1,zz),all_match(:,2,zz),1);
        all_match_sens(zz) = coeffs(1);

        %  store example data for plotting in figure
        if sess==1 & strfind(all_sess_files(zz).name,'6PG6')                
            summary.example.vismat_model = vismat;
            summary.example.vismat_data = Z.hexa_data_an.visits;

            vismat_t = find(sum(summary.example.vismat_data,1)==1);
            
            %   compute the distribution of cumulative visits and overlay
            for iii = 1:num_iter
                all_vis_cum(:,:,iii) = cumsum(summary.example.vismat_model(:,:,iii),2);
            end
            summary.example.cumvis.mean     = squeeze(mean(all_vis_cum,3));
            summary.example.cumvis.std      = squeeze(std(all_vis_cum,[],3));
            
            grima = TNC_CreateRBColormap(8,'grima');
            
            hhh = figure(500); clf;
            for jjj=1:6
                shadedErrorBar(vismat_t,summary.example.cumvis.mean(jjj,:),summary.example.cumvis.std(jjj,:),{'color',grima(jjj,:)});
                hold on;
            end

            for jjj=1:6
                plot(vismat_t,cumsum(summary.example.vismat_data(jjj,vismat_t)),'k-','Color',grima(jjj,:),'LineWidth',3);            
            end
            exportgraphics(hhh, '6PG6-ExampleCumVis.eps',"ContentType","vector");
        end

        % run for dopamine alpha fit
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(0.001,best_a2,best_a2,f.a4,f.a5,'sig_exp',S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,3) = mean(trans_r2_iter);
        opt.labels{3} = 'AQUA DA=alpha';
        opt.rew(zz,3) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,3) = mean(income_r2_iter);

        % run AQUA with flat alpha and no distance
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(com_in_param_space(1),0,0,com_in_param_space(2),com_in_param_space(3),'sig_exp',S.visit_matrix,ones(size(S.cost_per_port)),S.rew_sched,S.income,S.prior);
        end

        opt.r2(zz,4) = mean(trans_r2_iter);
        opt.labels{4} = 'AQUA -dist -dyn';
        opt.rew(zz,4) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,4) = mean(income_r2_iter);
        
        % run AQUA without distance weighting
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(0.001,com_in_param_space(1),com_in_param_space(1),com_in_param_space(2),com_in_param_space(3),'sig_exp',S.visit_matrix,ones(size(S.cost_per_port)),S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,5) = mean(trans_r2_iter);
        opt.labels{5} = 'AQUA -dist';
        opt.rew(zz,5) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,5) = mean(income_r2_iter);

        % run Q learn static alpha optimum
        alpha = [mean(alpha_vis) mean(alpha_vis)]; 
        beta = 1; % Softmax term 
        % Grossman model uses ~ Beta = 4, alpha(1) = 0.5, alpha(2) = 0.1
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_Q(alpha,beta,S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,6) = mean(trans_r2_iter);
        opt.labels{6} = 'Q-learn';
        opt.q.params.alpha = alpha; 
        opt.q.params.beta = beta; 
        opt.rew(zz,6) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,6) = mean(income_r2_iter);

        % run WSLS
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_WSLS(S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,7) = mean(trans_r2_iter);
        opt.labels{7} = 'WSLS';
        opt.rew(zz,7) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,7) = mean(income_r2_iter);
        
        % run OIO
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_OIO(S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,8) = mean(trans_r2_iter);
        opt.labels{8} = 'OIO';
        opt.rew(zz,8) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,8) = mean(income_r2_iter);
        
        % run RANDOM
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_RAND(S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,9) = mean(trans_r2_iter);
        opt.labels{9} = 'Random';
        opt.rew(zz,9) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,9) = mean(income_r2_iter);

% ------------- 
% ------------- RUN RANGE OF MODEL SIMULATIONS

        all_recloc  = [all_recloc numel(strfind(all_sess_files(zz).name,'NAc'))];
        all_taus    = [all_taus f.a5];

        if numel(strfind(all_sess_files(zz).name,'NAc'))

            figure(9);
            subplot(5,4,zz);

            plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',targety,'.','color', [0.8 0.8 0.8] );
            hold on;
            plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',movmean(targety,21), 'k-')
            plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',f([1:numel(Z.hexa_data_an.da_resp_all.r)]'),'r','linewidth',3);
            
            % plot(f,x(round(end/20):end),y(round(end/20):end)); title(num2str(-1/f.d)); 
            axis([0 500 -0.25 1.25]);

            figure(11);        
            subplot(5,4,zz);
            da_rew_resp = Z.hexa_data_an.da_resp_data.wins(Z.hexa_data_an.da_resp_data.r_vis_id==1,:);
            plot(Z.hexa_data_an.da_resp_data.range,mean(da_rew_resp(1:round(end/3),:))); hold on;
            plot(Z.hexa_data_an.da_resp_data.range,mean(da_rew_resp(round(end/3):round(2*end/3),:))); hold on;
            plot(Z.hexa_data_an.da_resp_data.range,mean(da_rew_resp(round(2*end/3):end,:))); hold on;
            plot(Z.hexa_data_an.da_resp_data.range,mean(Z.hexa_data_an.da_resp_data.wins(Z.hexa_data_an.da_resp_data.r_vis_id==0,:)));
            axis([min(Z.hexa_data_an.da_resp_data.range) max(Z.hexa_data_an.da_resp_data.range) -1 4]); box off;
            title(all_sess_files(zz).name(1:10));
            drawnow;

        end

        close(Z.hexa_data_an.da_hand1);
        close(Z.hexa_data_an.da_hand2);

end


session(sess).opt       = opt;

session(sess).all_coms  = all_coms;
session(sess).all_taus  = all_taus;
session(sess).all_r2    = all_r2;
session(sess).all_recloc= all_recloc;

%% Look at boxchart for the new runs

figure(57); clf;
boxchart(opt.r2(all_recloc==1,[1:2 6:9]).^2);
xticklabels(opt.labels([1:2 6:9]));
ylabel('Transition matrix similarity (r^2)');
xlabel('Model type');
ylim([-0.1 1]);

[p,t,stats] = anova1(opt.r2(all_recloc==1,[1:2 6:9]));
[c,m,h,~] = multcompare(stats);

figure(58); clf;
boxchart(opt.rew(all_recloc==1,[1:2 6:9])-opt.rew_act(all_recloc==1));
xticklabels(opt.labels([1:2 6:9]));
ylabel('\Delta Predicted Rewards Collected');
xlabel('Model type');
ylim([-200 100]);

[p_rew,t_rew,stats_rew] = anova1(opt.rew(:,[1:2 6:9])-opt.rew_act);
[c_rew,m_rew,h_rew,~] = multcompare(stats_rew);

figure(59); clf;
subplot(121);
boxchart(opt.r2(all_recloc==1,2:3).^2);
xticklabels(opt.labels([2:3]));
ylabel('Transition matrix similarity (r^2)');
ylim([0 1]);

subplot(122);
boxchart(opt.rew(all_recloc==1,[2:3])-opt.rew_act(all_recloc==1));
xticklabels(opt.labels([2:3]));
ylabel('\Delta Predicted Rewards Collected');
ylim([-200 100]);


figure(580); clf;
boxchart(opt.r2(:,[2 5 4]).^2);
xticklabels(opt.labels([2 5 4]));
ylabel('Transition matrix similarity (r^2)');
xlabel('Model type');
ylim([-0.1 1]);

[p,t,stats_aq] = anova1(opt.r2(:,[2 5 4]).^2);
[c_aqua,m_aqua,h_aqua,~] = multcompare(stats_aq)

sessions(sess).c_r2 = c;
sessions(sess).c_rew = c_rew;
sessions(sess).c_aqua = c_aqua;

%% Figure panel plotting routine for example session Data and Model transiiton matrices and cumulative visits

% General idea, but vismat is a little different structure
% tmp                     = find(sum(visit_matrix,1)==1);
% [~,visit_list_data]     = max(visit_matrix(:,tmp),[],1);    
% [trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data(1:end),0,1);
% [~,visit_list_model]    = max(hexa_model.visits(:,tmp),[],1);
% [trans_mat_model]       = HX_ComputeTransitionMatrix(visit_list_model(1:end),0,1);

vismat_t = find(sum(summary.example.vismat_data,1)==1);

%   compute the distribution of cumulative visits and overlay
for iii = 1:num_iter
    all_vis_cum(:,:,iii) = cumsum(summary.example.vismat_model(:,:,iii),2);
end
summary.example.cumvis.mean     = squeeze(mean(all_vis_cum,3));
summary.example.cumvis.std      = squeeze(std(all_vis_cum,[],3));

grima = TNC_CreateRBColormap(8,'grima');

figure(500); clf;
for jjj=1:6

    shadedErrorBar(vismat_t,summary.example.cumvis.mean(jjj,:),summary.example.cumvis.std(jjj,:),{'color',grima(jjj,:)});
    hold on;
    plot(vismat_t,cumsum(summary.example.vismat_data(jjj,vismat_t)),'k-','Color',grima(jjj,:),'LineWidth',3);

end

exag = TNC_CreateRBColormap(100,'exag');

[~,visit_list_data]     = max(summary.example.vismat_data(:,vismat_t),[],1);    
[trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data(1:end),0,1);

for iii = 1:num_iter
    [~,visit_list_model]    = max(summary.example.vismat_model(:,:,iii),[],1);
    trans_mat_model(:,:,iii)= HX_ComputeTransitionMatrix(visit_list_model(1:end),0,1);
end

figure(501); clf;
subplot(121);
imagesc(trans_mat_data,[0 0.2]); colormap(exag);
axis equal; axis off; title('Data');
subplot(122);
imagesc(squeeze(mean(trans_mat_model,3)),[0 0.2]); colormap(exag);
axis equal; axis off; title('Model');


%% Figure panel for compring (r2/nLL) and matching sensitivity from session 1 to session 2

[hexa_map] = TNC_CreateRBColormap(6,'grima');    
figure(19);

subplot(1,2,sess);

for zz=1:size(all_match,3)
    if zz==1
        hold off;
        plot([-3 0],[-3 0],'k--');
        hold on;
    end
    scatter(all_match(:,1,zz),all_match(:,2,zz),100,1:6,'filled'); colormap(hexa_map);    
    axis([-2.25 0 -2.25 0]);
end

title(['Session ' num2str(sess) ' : s= ' num2str(mean(all_match_sens))]);

session(sess).all_match = all_match;
session(sess).all_match_sens = all_match_sens;

%% Sensitivity subpanel

figure(61);  clf;
boxchart(1+[zeros(size(session(1).all_match_sens)) ; ones(size(session(2).all_match_sens))],[session(1).all_match_sens ; session(2).all_match_sens]);
hold on;
scatter(1.4+[zeros(size(session(1).all_match_sens)) ; ones(size(session(2).all_match_sens))],[session(1).all_match_sens ; session(2).all_match_sens],50,'k');
axis([0.5 2.5 0 1]);
ylabel('Sensitivity'); xlabel('Session');

p = anova1([session(1).all_match_sens ; session(2).all_match_sens],[zeros(size(session(1).all_match_sens)) ; ones(size(session(2).all_match_sens))]);

figure(61);
text(1.5,0.95,['p=' num2str(p)]);

%% Sensitivity subpanel
figure(65); clf;
boxchart(1+[zeros(size(session(1).opt.r2(:,2))) ; ones(size((session(2).opt.r2(:,2))))],[session(1).opt.r2(:,2) ; session(2).opt.r2(:,2)].^2);
axis([0.5 2.5 0 1]);
ylabel('Transition matrix similarity (r^2)');
xlabel('Session');

p = anova1([session(1).opt.r2(:,2).^2 ; session(2).opt.r2(:,2).^2],[zeros(size(session(1).opt.r2(:,2))) ; ones(size((session(2).opt.r2(:,2))))]);

figure(65);
text(1.5,0.95,['p=' num2str(p)]);


%% Plotting routines for examining model comparisons
% 
% figure(900); clf;
% plot([250 750],[250 750],'k-'); hold on;
% scatter(all_coms(all_recloc==1,3),all_taus(all_recloc==1),100,all_r2(all_recloc==1,1)); colormap(sym);
% scatter(mean(all_coms(all_recloc==1,3)),mean(all_taus(all_recloc==1)),100,'k','filled');
% axis([100 750 100 750]);

figure(850); clf;
figure(851); clf;

com_labels = {'Magnitude','Rise','Decay'};
for sess=1:2

    figure(850);
    for bb=1:3
        subplot(1,3,bb);
        hold on;
        boxchart(sess.*ones(size(session(sess).all_coms(:,bb))),session(sess).all_coms(:,bb)); hold on;

        if sess==1
            title(com_labels{bb});
        end

        [p_val] = ranksum(session(1).all_coms(:,bb),session(2).all_coms(:,bb));
        title(['P-value: ' num2str(p_val)]);
    end

    figure(851);
    boxchart(sess*ones(size(session(sess).all_taus(session(sess).all_recloc==1))),session(sess).all_taus(session(sess).all_recloc==1)); hold on;
    
end



%% New optimization script that also allows for optimization of distance scaling

clear model_compare


all_files = dir('~/Dropbox (HHMI)/hexaport/photometry/full_dataset/*conc_b*');
% all_files = dir('~/Dropbox (HHMI)/hexaport/photometry/full_dataset/6PG31_NAc_conc_b*');

path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/photometry/full_dataset/';


pathcost_logic  = 1
session         = 5
notes = ['da_store_analyzed_sess' num2str(session) 'fitAlpha_fitPathCost'];
dir_path = [notes '/']
[SUCCESS,~,~] = mkdir(path,dir_path);

cost_per_port =                 ...
[1	14	18	70	72.2	65.5;   ...
14	1	22.8	56	65.5	42; ...
18	22.8	1	72.2	70	56; ...
70	56	72.2	1	18	22.8;   ...
72.2	65.5	70	18	1	14; ...
65.5	42	56	22.8	14	1];

if pathcost_logic
    % pathcost = cost_per_port;
    % pathcost(pathcost==14)=7;
    % pathcost(pathcost==18)=9;
    % pathcost(pathcost==56)=28;
    % pathcost(pathcost==70)=35;

    % if shelter acts like a wall
    pathcost = cost_per_port./2;
    pathcost(pathcost==11.4)=22.8; % c is the only path that doesn't go along a wall nor through a shelter
    pathcost(pathcost==0.5)=1; % c is the only path that doesn't go along a wall nor through a shelter

    cost_per_port=pathcost
end

photo_flag = 1;
figure(600);

port_color_map      = TNC_CreateRBColormap(8,'mapb');

testing = 1;
for mmm = 1:numel(all_files) % mice 11 and 16 do not have session 2 data
    
    clear prior;
    breaks = strfind(all_files(mmm).name,'_');
    mouse_name = all_files(mmm).name(1:breaks(1)-1);
      
    photo_filename  = [path mouse_name '_photo.csv'];
    [hexa_data]     = HX_load_csv([path all_files(mmm).name], 0, photo_flag, photo_filename);

    if min(ismember(session,unique(hexa_data.session_n)))==1

            [hexa_data_an]  = HX_analyze_session(hexa_data,session,photo_flag);
        
            intervals = [30 60 240 1200 2400];
            port_intervals = zeros(numel(session),6);
            for ss=session
                for qq=1:6
                    port_intervals(ss,qq)       = intervals(unique(hexa_data.port_rank(hexa_data.port_n==qq & ismember(hexa_data.session_n,ss))));
                    port_rank_this_sess(ss,qq)  = (unique(hexa_data.port_rank(hexa_data.port_n==qq & ismember(hexa_data.session_n,ss))));
                end
            end    
        
        %----------------------------------------
        %----------------------------------------
        %-------------- FITTING ALPHA PER SESSION
        %----------------------------------------
        %----------------------------------------

        visit_matrix = hexa_data_an.visits;

        sample_logic = sum(visit_matrix,1);
        
        all_visits = find(sum(hexa_data_an.visits,1)==1);
        rew_logic = sum(hexa_data_an.rewards,1);
        all_rewards = rew_logic(all_visits);
        income = movmean(all_rewards,51);
        
        % sampling rate is now set to 1 Hz
        frame_rate = 1;
        
        hexa_model.rew_sched = zeros(size(hexa_data_an.visits));
        for ss=unique(hexa_data_an.sessID)'
            valid_inds = find(hexa_data_an.sessID==ss);
            for qq=1:6
                hexa_model.rew_sched(qq,valid_inds(1):round(port_intervals(ss,qq)*frame_rate):valid_inds(end)) = 1;
            end
        end
        hexa_model.rew_sched(:,2) = 1;
        
        rew_sched = hexa_model.rew_sched;
        
        % example values for alpha_params
        alpha_params_init = [0 0.2 0.2 300 500];
        v_ind = 1:sum(sample_logic);

        vis_inds                = find(sum(visit_matrix,1)==1);
        [~,visit_list_data]     = max(visit_matrix(:,vis_inds),[],1);    
        [trans_mat_data]        = HX_ComputeTransitionMatrix(visit_list_data(1:end),0,1);
        
        % Visualize the alpha function with default
        % alpha_version = 'doub_exp'
        % alpha_vis_init = alpha_params_init(1) + (alpha_params_init(2)*(1-exp(-v_ind/alpha_params_init(4))) .* (alpha_params_init(3)*exp(-v_ind/alpha_params_init(5))));
        
        alpha_version = 'sig_exp'        
        [rise_kern] = TNC_CreateGaussian(alpha_params_init(4),100,sum(sample_logic),1);
        alpha_rise = cumsum(rise_kern)*alpha_params_init(2);
        alpha_vis_init = alpha_params_init(1) + (alpha_rise .* (alpha_params_init(3)*exp(-v_ind/alpha_params_init(5))));
        
        % Online plot initialization
        exag_map    = TNC_CreateRBColormap(8,'exag');
        figure(10); clf; 
        figure(11); clf;
        
        % Examine the base r2 with a guess at alpha initialization
        a1 = alpha_params_init(1);
        a2 = alpha_params_init(2);
        a3 = alpha_params_init(3);
        a4 = alpha_params_init(4);
        a5 = alpha_params_init(5);
        
        num_iter        = 8;
        vismat          = zeros(6,numel(all_visits),num_iter);
        rewmat          = zeros(6,numel(all_visits),num_iter);
        trans_mat_modrun= zeros(6,6,num_iter);

        if session==1
            % assume random initialization
            prior(:,1) = 1/6;
            prior(:,2) = 0.1;
        else
            % estimate prior from hexa_data_an for session-1
            [hexa_data_an_prior]    = HX_analyze_session(hexa_data,session-1,0);
            tmp                     = find(sum(hexa_data_an_prior.visits,1)==1);
            [~,visit_list_data]     = max(hexa_data_an_prior.visits(:,tmp(round(0.7*end:end))),[],1);        
            [trans_mat_data_prior]  = HX_ComputeTransitionMatrix(visit_list_data,0,1);
            prior(:,2)              = diag(trans_mat_data_prior,0);
            prior(:,1)              = sum(trans_mat_data_prior,1)' - diag(trans_mat_data_prior,0);
            % add noise?
        end


        %----------------------------------------
        %-------------- FITTING ALPHA TO DA
        %----------------------------------------
        alpha = @(a1,a2,a3,a4,a5,x) a1 + (a2 ./ (1+exp((a4-x)/(a4./6)))) .*  (a3*exp(-x/a5));
        fitfun = fittype( alpha );        
        targety = movmean(hexa_data_an.da_resp_all.r,3)./max(movmean(hexa_data_an.da_resp_all.r,11));
        a0 = [ 0 0.5 0.5 100 1000 ];        
        [dopa.f,dopa.gof] = fit([1:numel(hexa_data_an.da_resp_all.r)]',targety,fitfun,'StartPoint',a0,'Upper',[0.1 1 1 numel(targety) 2*numel(targety)],'Lower',[0 0 0 20 20]);

        %----------------------------------------
        %-------------- FITTING ALPHA TO DA
        %----------------------------------------
        % grid search optimization
        a2_vec = [0.05 0.1 0.2 0.5 0.99];
        a4_vec = [10 20 50 100 200];
        a5_vec = [50 100 200 500 1000];
        de_vec = [0.5 1 1.5 2 3];

        [~,close_a4] = min(abs(dopa.f.a4-a4_vec));
        [~,close_a5] = min(abs(dopa.f.a5-a5_vec));

        for a2 = a2_vec
            for a4 = a4_vec
                % for a5 = a5_vec & original
                for de = de_vec
        
                    a1 = 0.001;
                    a2 = a2;
                    a3 = a2;

                    a5 = 5 * a4; % some function of a4

                    trans_r2_iter  = zeros(1,num_iter);
                    income_r2_iter = zeros(1,num_iter);

                    parfor iter = 1:num_iter
                        [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(a1,a2,a3,a4,a5,alpha_version,visit_matrix,cost_per_port.^de,rew_sched,income,prior);
                    end
        
                    opt_r2_tensor(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = median(trans_r2_iter);
                    opt_inc_tensor(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = median(income_r2_iter);
                    params_a2(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = a2;
                    params_a4(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = a4;
                    params_de(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = de;

                    LL = sum( (reshape(trans_mat_data+0.001,1,36)) .* log( reshape(trans_mat_data+0.001,1,36) ./ reshape(squeeze(mean(trans_mat_modrun,3)+0.001),1,36) ) );
                    opt_LL_tensor(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = mean(LL);

                    tot_rew     = sum(sum(mean(rewmat,3)));               
                    opt_RColl_tensor(find(a2==a2_vec),find(a4==a4_vec),find(de==de_vec)) = tot_rew;
                    
                    if find(de==de_vec)==numel(de_vec)
                        hhh = figure(11); 
                        subplot(4,numel(a2_vec),find(a2==a2_vec));
                        imagesc(squeeze(opt_r2_tensor(find(a2==a2_vec),:,:)),[0.25 0.85]); colormap(exag_map);                        
                        title('Trans R2');
                        figure(11); hold off; subplot(4,numel(a2_vec),find(a2==a2_vec)+numel(a2_vec));
                        imagesc(squeeze(opt_inc_tensor(find(a2==a2_vec),:,:)),[0.05 0.25]); colormap(exag_map);
                        hold on;
                        scatter(close_a5,close_a4,50,'k','filled');
                        title('Income RMSE');
                        figure(11); subplot(4,numel(a2_vec),find(a2==a2_vec)+numel(a2_vec)+numel(a2_vec));
                        imagesc(squeeze(opt_LL_tensor(find(a2==a2_vec),:,:)),[0.5 5]); colormap(exag_map);
                        title('nLL');
                        figure(11); subplot(4,numel(a2_vec),find(a2==a2_vec)+numel(a2_vec)+numel(a2_vec)+numel(a2_vec));
                        imagesc(abs(squeeze(opt_RColl_tensor(find(a2==a2_vec),:,:))-sum(sum(hexa_data_an.rewards))),[0 500]); colormap(exag_map);
                        title(['Total Reward - mouse: ' num2str(sum(sum(hexa_data_an.rewards)))]);
                    end
                end
            end
        end

        exportgraphics(hhh, [path dir_path mouse_name '_FitTensSummary.pdf'],"ContentType","vector"); % write out Fig 11
        
        % look for joint min
        summary_fit_fig = figure(700);
        subplot(ceil(numel(all_files)/5),5,mmm);
        scatter(reshape(opt_inc_tensor,1,5*5*5),reshape(opt_r2_tensor,1,5*5*5),(1+reshape(params_a2,1,5*5*5)).^2*50,reshape(params_de,1,5*5*5),'filled','MarkerEdgeColor','k'); colormap(exag_map);
        xlabel('inc RMSE'); ylabel('Transition matrix r^2');
        title([mouse_name '; r2: ' num2str(max(reshape(opt_r2_tensor,1,5*5*5))) '; rmse: ' num2str(min(reshape(opt_LL_tensor,1,5*5*5)))]);
        axis([0.05 0.25 0 1]);
        
        %----------------------------------------
        %----------------------------------------
        %-------------- FITTING ALPHA PER SESSION
        %----------------------------------------
        %----------------------------------------
        
           % interim visualization
            Nback               = 1;
            hh = figure(499); clf;
            subplot(1,3,1:2);
            % plot(hexa_data.e√vent_time_con,[0 diff(hexa_data.session_n)'],'color',[0 0 0 0.25],'linewidth',2); hold on;
            for pp=1:6
                plot(hexa_data.event_time_con(hexa_data_an.visit_indices),hexa_data_an.p_choice_all(pp,:),'linewidth',2,'color',port_color_map(unique(hexa_data.port_rank(hexa_data.port_n==pp & ismember(hexa_data.session_n,ss))),:)); hold on;
            end
            axis([0 max(hexa_data.event_time_con(hexa_data_an.visit_indices)) 0 0.67]);
            box off;
            ylabel('P(visit,port)'); 
            title(mouse_name);
            
            [~,ii_r] = sort(port_rank_this_sess(1,:));
        
            tmp                 = find(sum(hexa_data_an.visits,1)==1);
            [~,visit_list_data] = max(hexa_data_an.visits(ii_r,tmp),[],1);    
        
            [trans_mat_data]    = HX_ComputeTransitionMatrix(visit_list_data,0,Nback);
            exag = TNC_CreateRBColormap(8,'exag');
        
            subplot(133);
            imagesc(trans_mat_data,[0 0.2]); 
            colormap(exag); axis equal; box off; colorbar;
            title([mouse_name '; Nback=' num2str(Nback)]);
            drawnow;
        
            summary_fig = figure(600);    
            subplot(ceil(numel(all_files)/5),5,mmm)
            imagesc(trans_mat_data,[0 0.25]); colormap(exag); axis equal; box off; colorbar;
            title([ mouse_name '; r# ' num2str(sum(sum(hexa_data_an.rewards,1))) '; v# ' num2str(sum(sum(hexa_data_an.visits,1))) '; r2: ' num2str(max(reshape(opt_r2_tensor,1,5*5*5)))]);
            drawnow;
        
            eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_sess' num2str(session) '_opt.mat a* de* opt* dopa tot_rew visit_matrix cost_per_port rew_sched income prior port_rank_this_sess']);
            eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_sess' num2str(session) '_an.mat hexa_data_an']);
            disp(['Completed fitting for ' all_files(mmm).name ' session(s): ' num2str(session)]);

    else

        disp(['Skipped ' mouse_name ' because data lacked one of the target sessions; ' 'target session array: ' num2str(session) ' sessions found: ' num2str(unique(hexa_data.session_n)')])

    end

end

exportgraphics(summary_fig, [path dir_path 'all_mouse_summary_trans_mat.pdf'],"ContentType","vector");
exportgraphics(summary_fit_fig, [path dir_path 'all_mouse_summary_fit_mat.pdf'],"ContentType","vector");

%% Take saved files and calculate summary stats for session 2-3 transition analysis

a=1;
b=0;
c=0;
frac = 0.95;

figure(9); clf;
figure(11); clf;
figure(12); clf;
figure(13); clf;

for sess=1:3

figure(799+sess); clf;    
    clear opt

all_coms    = [];
all_isos    = [];
all_taus    = [];
all_a2s     = [];
all_recloc  = [];
all_r2      = [];
all_nLL     = [];
all_mrew    = [];

all_sess_files = dir(['*sess' num2str(sess) '*_opt.mat']);
opt.r2      = zeros(numel(all_sess_files),9); % alpha fit, com fit, da fit, da raw, q, oio, wsls, random
opt.inc     = zeros(numel(all_sess_files),9); % alpha fit, com fit, da fit, da raw, q, oio, wsls, random
opt.rank    = zeros(numel(all_sess_files),1);
opt.rew     = zeros(numel(all_sess_files),1);

all_match       = zeros(6,2,numel(all_sess_files));
all_match_sens  = zeros(numel(all_sess_files),1);

for zz=1:numel(all_sess_files)

    clear S Z
    S = load(all_sess_files(zz).name);

    breaks = strfind(all_sess_files(zz).name,'_');

    loss = a*S.opt_r2_tensor-b*S.opt_inc_tensor-c*S.opt_LL_tensor;
    
    [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));
    sym = TNC_CreateRBColormap(numel(top_xperc_inds),'cpb');
        
    [a2_inds,a4_inds,de_inds] = ind2sub(size(S.opt_r2_tensor),top_xperc_inds);
    
% forgot to save, kludge    
S.de_vec = [0.5 1 1.5 2 3];

    [com_in_param_space] = centerOfMass3D(S.a2_vec(a2_inds), S.a4_vec(a4_inds), S.de_vec(de_inds), loss(top_xperc_inds)');

    all_r2      = [all_r2 ; max(S.opt_r2_tensor,[],"all")];
    all_coms    = [all_coms ; com_in_param_space];
    all_isos    = [all_isos ; numel(top_xperc_inds)];

    % load the corresponding hexa_data_an struct
    Z = load([all_sess_files(zz).name(1:end-7) 'an.mat']);

    % Get data/optimal discovered fit observations for comparison
    opt.r2(zz,1) = max(S.opt_r2_tensor,[],"all");
    opt.labels{1} = 'AQUA opt grid';
    opt.rew(zz,1) = S.tot_rew;
    opt.rew_act(zz,1) = numel(find(Z.hexa_data_an.rewards==1));

% ------------- ALPHA fitting routine

        alpha = @(a1,a2,a3,a4,a5,x) a1 + (a2 ./ (1+exp((a4-x)/(a4./6)))) .*  (a3*exp(-x/a5));
        fitfun = fittype( alpha );

        targety = movmean(Z.hexa_data_an.da_resp_all.r,3)./max(movmean(Z.hexa_data_an.da_resp_all.r,11));

        % reasonable initial guesses
        a0 = [ 0 0.5 0.5 100 1000 ];

        [f,gof] = fit([1:numel(Z.hexa_data_an.da_resp_all.r)]',targety,fitfun,'StartPoint',a0,'Upper',[0.1 1 1 numel(targety) 2*numel(targety)],'Lower',[0 0 0 10 20]);        

% ------------- ALPHA fitting routine

if sess==3 & numel(strfind(all_sess_files(zz).name,'ML15'))
        all_recloc  = [all_recloc 0];
else
        all_recloc  = [all_recloc numel(strfind(all_sess_files(zz).name,'NAc'))];
end

        all_taus    = [all_taus f.a5];
        all_a2s     = [all_a2s f.a2];
        all_mrew    = [all_mrew median(Z.hexa_data_an.da_resp_all.r)];

        if numel(strfind(all_sess_files(zz).name,'NAc'))

            figure(9);
            subplot(5,4,zz);

            plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',targety,'.','color', [0.8 0.8 0.8] );
            hold on;
            plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',movmean(targety,21), 'k-')
            plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',f([1:numel(Z.hexa_data_an.da_resp_all.r)]'),'r','linewidth',3);
            
            % plot(f,x(round(end/20):end),y(round(end/20):end)); title(num2str(-1/f.d)); 
            axis([0 500 -0.25 1.25]);

            figure(10+sess);        
            subplot(5,4,zz);
            da_rew_resp = Z.hexa_data_an.da_resp_data.wins(Z.hexa_data_an.da_resp_data.r_vis_id==1,:);
            plot(Z.hexa_data_an.da_resp_data.range,mean(da_rew_resp(1:round(end/3),:))); hold on;
            plot(Z.hexa_data_an.da_resp_data.range,mean(da_rew_resp(round(end/3):round(2*end/3),:))); hold on;
            plot(Z.hexa_data_an.da_resp_data.range,mean(da_rew_resp(round(2*end/3):end,:))); hold on;
            plot(Z.hexa_data_an.da_resp_data.range,mean(Z.hexa_data_an.da_resp_data.wins(Z.hexa_data_an.da_resp_data.r_vis_id==0,:)));
            axis([min(Z.hexa_data_an.da_resp_data.range) max(Z.hexa_data_an.da_resp_data.range) -1 4]); box off;
            title(all_sess_files(zz).name(1:10));
            drawnow;

        end

        close(Z.hexa_data_an.da_hand1);
        close(Z.hexa_data_an.da_hand2);

        session(sess).file(zz).name = all_sess_files(zz).name;

end


session(sess).opt       = opt;

session(sess).all_coms  = all_coms;
session(sess).all_taus  = all_taus;
session(sess).all_a2s   = all_a2s;
session(sess).all_mrew  = all_mrew;
session(sess).all_r2    = all_r2;
session(sess).all_recloc= all_recloc;
end

%% Plotting routines for examining model comparisons

figure(850); clf;
figure(851); clf;

com_labels = {'\alpha Mag.','\alpha tau','Path Exp.'};
for sess=1:3

    figure(850);
    for bb=1:3
        subplot(1,3,bb);
        hold on;
        boxchart(sess.*ones(size(session(sess).all_coms(:,bb))),session(sess).all_coms(:,bb)); hold on;

        if sess==1
            title(com_labels{bb});
        end

        [p_val] = ranksum(session(1).all_coms(:,bb),session(2).all_coms(:,bb));
        [p_val3] = ranksum(session(2).all_coms(:,bb),session(3).all_coms(:,bb));
        title(['P-value: ' num2str(p_val) ' P3-value: ' num2str(p_val3)]);
        ylabel(com_labels{bb}); xlabel('session');
    end

    figure(851);
    subplot(131);
    boxchart(sess*ones(size(session(sess).all_coms(session(sess).all_recloc==1,1))),session(sess).all_coms(session(sess).all_recloc==1,1)); hold on;

    title(['\alpha Mag.; p=' num2str(p_val3)]); xlabel('Session');

    subplot(132);
    boxchart(sess*ones(size(session(sess).all_mrew(session(sess).all_recloc==1))),session(sess).all_mrew(session(sess).all_recloc==1)); hold on;
    [p_mrew] = ranksum(session(2).all_mrew(session(2).all_recloc==1),session(3).all_mrew(session(3).all_recloc==1));
    title(['Med. Rew. Resp. ; p=' num2str(p_mrew)]);
    
end

%% Look at paired DA

for sess=2:3
    disp(sess);
    inds = find(session(sess).all_recloc==1);
    for jj=inds
        disp([num2str(jj) ': ' session(sess).file(jj).name]);
    end
end

update_sess2_recloc = zeros(1,16);
update_sess2_recloc([3 7 9 10 13 14 16])=1

    figure(851);
    subplot(133);
    boxchart([session(2).all_mrew(update_sess2_recloc==1)',session(3).all_mrew(session(3).all_recloc==1)']); hold on;

    title(['Rew. Resp. Paired Animals; p=' num2str(p_da)]);

[h_da,p_da] = ttest(session(2).all_mrew(update_sess2_recloc==1),session(3).all_mrew(session(3).all_recloc==1))