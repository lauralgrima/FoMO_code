% HX_OptoSim_Script
%% SCRIPT TO FIT OPTIMAL ALPHA FUNCTION BEFORE TESTING OPTOSTIM PREDICTIONS

clear model_compare

all_files = dir('~/Dropbox (HHMI)/hexaport/optogenetics/*4worst*');
path = '/Users/dudmanj/Dropbox (HHMI)/hexaport/optogenetics/';

cost_per_port =                 ...
[1	14	18	70	72.2	65.5;   ...
14	1	22.8	56	65.5	42; ...
18	22.8	1	72.2	70	56; ...
70	56	72.2	1	18	22.8;   ...
72.2	65.5	70	18	1	14; ...
65.5	42	56	22.8	14	1];

for session         = 1:2

    clear model_compare

    notes = ['da_store_analyzed_sess' num2str(session) 'nLL_OptoSched_Final'];
    dir_path = [notes '/']
    [SUCCESS,~,~] = mkdir(path,dir_path);
    
    photo_flag = 0;
    figure(600);
    
    port_color_map      = TNC_CreateRBColormap(8,'mapb');
    
    testing = 1;
    for mmm = 1:numel(all_files) % mice 11 and 16 do not have session 2 data
        
        clear prior;
        breaks = strfind(all_files(mmm).name,'_');
        mouse_name = all_files(mmm).name(1:breaks(1)-1);
        [hexa_data]     = HX_load_csv([path all_files(mmm).name], 0, photo_flag, '');
    
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
            reward_matrix = hexa_data_an.rewards;

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
            
            alpha_version = 'single_exp';
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
            
            num_iter        = 10;
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
            % grid search optimization
            a1_vec = [0.025 0.05 0.1 0.2 0.5];
            a2_vec = [0.025 0.05 0.1 0.2 0.5];
            a5_vec = [10 50 100 200 500];

            close_a4 = [];
            close_a5 = [];

            for a1 = a1_vec
                for a2 = a2_vec
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
                        opt_r2_tensor(find(a1==a1_vec),find(a2==a2_vec),find(a5==a5_vec)) = median(trans_r2_iter);
                        opt_inc_tensor(find(a1==a1_vec),find(a2==a2_vec),find(a5==a5_vec)) = median(income_r2_iter);

                        params_a1(find(a1==a1_vec),find(a2==a2_vec),find(a5==a5_vec)) = a1;
                        params_a2(find(a1==a1_vec),find(a2==a2_vec),find(a5==a5_vec)) = a2;
                        params_a5(find(a1==a1_vec),find(a2==a2_vec),find(a5==a5_vec)) = a5;
    
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
    
            % exportgraphics(hhh, [path dir_path mouse_name '_FitTensSummary.pdf'],"ContentType","vector"); % write out Fig 11
            
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
                        
                summary_fig = figure(600);    
                subplot(ceil(numel(all_files)/5),5,mmm)
                imagesc(trans_mat_data,[0 0.25]); colormap(exag); axis equal; box off; colorbar;
                title([ mouse_name '; r# ' num2str(sum(sum(hexa_data_an.rewards,1))) '; v# ' num2str(sum(sum(hexa_data_an.visits,1))) '; r2: ' num2str(max(reshape(opt_r2_tensor,1,5*5*5)))]);
                drawnow;
            
                if photo_flag==0
                    dopa = [];
                end

                eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_sess' num2str(session) '_alphaOnly_opt.mat a* opt* dopa tot_rew visit_matrix reward_matrix cost_per_port rew_sched income prior port_rank_this_sess']);
                eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_sess' num2str(session) '_alphaOnly_an.mat hexa_data_an']);

                if session==1
                    eval(['save ~/Downloads/' all_files(mmm).name(1:end-4) '_AllSess' '_dat.mat hexa_data']);
                end
                
                disp(['Completed fitting for ' all_files(mmm).name ' session(s): ' num2str(session)]);
    
        else
    
            disp(['Skipped ' mouse_name ' because data lacked one of the target sessions; ' 'target session array: ' num2str(session) ' sessions found: ' num2str(unique(hexa_data.session_n)')])
    
        end
    
    end
    
    exportgraphics(summary_fig, [path dir_path 'all_mouse_summary_trans_mat.pdf'],"ContentType","vector");
    exportgraphics(summary_fit_fig, [path dir_path 'all_mouse_summary_fit_mat.pdf'],"ContentType","vector");
end


%% SCRIPT TO EXAMINE THE EFFECTS OF STIM ON INTERVAL SCHEDULE FOR FINAL FIGURE

%----------------------------------------
%----------------------------------------
% Let the data control which animal is being examined and then look over sessions for fit params within that    
%----------------------------------------
%----------------------------------------
all_sess_files = dir('*dat.mat');

figure(900); clf; figure(101); clf; 
unos = 1:3:34; dos = 2:3:35; tres = 3:3:36;
clear GLM_export mouse;

loss_target = 'combined'

for zz=1:numel(all_sess_files)

    %----------------------------------------
    %----------------------------------------
    % load the corresponding hexa_data_an     
    %----------------------------------------
    %----------------------------------------
    T = load(all_sess_files(zz).name);

    T.hexa_data.event_time_con(session_bounds:end) = T.hexa_data.event_time_con(session_bounds:end)+round(T.hexa_data.event_time_con(session_bounds-1));
    session_bounds      = find([0 diff(T.hexa_data.session_n')]==1);
    
    clear da_resp_rew da_resp_ure aqua_model all_com;

    %----------------------------------------
    %----------------------------------------
    % Examine the evolution of port choice over all sessions
    %----------------------------------------
    %----------------------------------------
    gmap                = TNC_CreateRBColormap(8,'grima');
    sess_map            = TNC_CreateRBColormap(6,'yb');
    sess_map            = sess_map(4:end,:);

    visit_matrix        = zeros(6,ceil(max(unique(T.hexa_data.event_time_con))));
    reward_matrix       = zeros(6,ceil(max(unique(T.hexa_data.event_time_con))));
    visit_matrix_sm     = zeros(6,ceil(max(unique(T.hexa_data.event_time_con))));
    reward_matrix_sm    = zeros(6,ceil(max(unique(T.hexa_data.event_time_con))));

    session_ids         = ones(1,ceil(max(unique(T.hexa_data.event_time_con))));
    choice_kernel       = [0 ones(1,10*60) 0]./(10*60);
    rew_kernel          = [0 ones(1,30*60) 0]./(30*60);

    if numel(session_bounds)==0
        sess_start_times    = [];
    else
        sess_start_times    = round(T.hexa_data.event_time_con(session_bounds))
    end
    
    % In this data the event times are not continuous for session 1 and 2
    % so I need to make them

    uv_inds             = find(T.hexa_data.unique_vis==1);
    times               = T.hexa_data.event_time_con(uv_inds);
    if min(diff(times))<1
        when = find([1 diff(times')]<1);
        while numel(when)>0
            when = find([1 diff(times')]<1);
            times(when) = times(when)+0.5;
            disp([num2str(numel(when)) ' shifts applied...']);
            when = find([1 diff(times')]<1);
        end
        disp(['Done... min(delta(time))=' num2str( min(diff(times)))]);
    end
    ports               = T.hexa_data.port_n(uv_inds);

    rw_inds             = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==1);
    rtimes              = T.hexa_data.event_time_con(rw_inds);
    rports              = T.hexa_data.port_n(rw_inds);

    %----------------------------------------
    %----------------------------------------
    % Get rank order of port quality per session
    %----------------------------------------
    %----------------------------------------
    per_sess_rank = zeros(6,numel(unique(T.hexa_data.session_n)));
    for jjj=unique(T.hexa_data.session_n)'
        for kkk=1:6
            uv_ps_inds = find(T.hexa_data.unique_vis==1 & T.hexa_data.session_n==jjj & T.hexa_data.port_n==kkk);
            per_sess_rank(kkk,jjj) = unique(T.hexa_data.port_rank(uv_ps_inds));
        end
    end
    
    if numel(sess_start_times)==0
        session_ids=ones(1,size(visit_matrix,2));
    else
        session = zeros(1,ceil(max(unique(T.hexa_data.event_time_con))));
        for jjj = 1:numel(sess_start_times)-1
            session_ids(sess_start_times(jjj):sess_start_times(jjj+1))=jjj+1;
        end
        session_ids(sess_start_times(end):end)=numel(sess_start_times)+1;
    end

    %----------------------------------------
    %----------------------------------------
    % use this per_sess_rank and session_ids and visit_matrix to make concatenated reward schedule for entire animal's data
    %----------------------------------------
    %----------------------------------------
    rew_sched = zeros(size(visit_matrix));
    intervals = [30 60 240 1200 2400];
    for jjj=unique(T.hexa_data.session_n)'
        valid_inds = find(session_ids==jjj);
        these_intervals = intervals(per_sess_rank(:,jjj));
        for qq=1:6
            rew_sched(qq,valid_inds(1):round(these_intervals(qq)):valid_inds(end)) = 1;
        end
    end
    rew_sched(:,2) = 1;
    
    figure(9); clf;
    for qqq = 1:6
        visit_matrix(qqq,ceil(times(ports==qqq)))=1;
        visit_matrix_sm(qqq,:) = movmean(visit_matrix(qqq,:),10*60);

        figure(9); subplot(211);
        plot(visit_matrix_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 size(visit_matrix,2)]); ylim([0 0.15]);
        ylabel('Fraction of visits');

        reward_matrix(qqq,ceil(rtimes(rports==qqq)))=1;
        reward_matrix_sm(qqq,:) = movmean(reward_matrix(qqq,:),30*60);

        figure(9); subplot(212);
        semilogy(reward_matrix_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 size(reward_matrix,2)]); ylim([0 0.075]);
        ylabel('Fraction of visits rewarded');

        figure(900);
        subplot(6,6,dos(zz));
        plot(visit_matrix_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 180*5*60]);box off;
        if zz==1
            ylabel('Fraction of visits'); 
        end
    end
    figure(9); subplot(211); yyaxis right; plot(session_ids,'k'); box off;
    all_rewards = sum(reward_matrix_sm,1);
    ylabel('Session');
    figure(9); subplot(212); plot(all_rewards,'k-'); yyaxis right; plot(session_ids,'k'); box off;
    ylabel('Session');

    all_rewards     = sum(reward_matrix,1);
    all_visits      = find(sum(visit_matrix,1)==1);
    all_sessid      = session_ids(all_visits);
    income          = movmean(all_rewards(all_visits),51);

    %----------------------------------------
    %----------------------------------------
    % Load all session data an find out the optimal fit params
    %----------------------------------------
    %----------------------------------------
    dir_breaks = strfind(T.hexa_data.filename,'/');
    just_file = T.hexa_data.filename(dir_breaks(end)+1:end);
    mname_end = strfind(just_file,'_');
    mouse_name = just_file(1:mname_end(1)-1);
    all_opt_files = dir([mouse_name '*opt*']);
    alpha = [];
    frac=0.98;

    all_com = zeros(numel(all_opt_files),3);
    alpha_da = @(a1,a2,a5,x) a1 + (a2*exp(-x/a5));

    for zzz=1:numel(all_opt_files)

        S = load(all_opt_files(zzz).name);
            
        switch loss_target

            case 'combined'
                % Loss function combines r2 and income RMSE
                loss = S.opt_r2_tensor - S.opt_inc_tensor;                
                [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));

            case 'trans'
                % Loss function combines r2 and income RMSE
                loss = S.opt_r2_tensor;                
                [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));

            case 'income'
                % Loss function jsut uses income RMSE
                loss = S.opt_inc_tensor;                
                [top_xperc_inds] = find(loss<(1+(1-frac))*min(loss,[],"all"));

            case 'freeze_de'
                % Loss function jsut uses income RMSE
                loss = S.opt_r2_tensor - S.opt_inc_tensor;             
                [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));
                
        end

        [a2_inds,a5_inds,a1_inds] = ind2sub(size(S.opt_r2_tensor),top_xperc_inds);
        [all_com(zzz,:)] = centerOfMass3D(S.a2_vec(a2_inds), S.a5_vec(a5_inds), S.a1_vec(a1_inds), loss(top_xperc_inds)');

        v_ind = 1:sum(sum(S.visit_matrix,1));
        
        if zzz==1
            % use sig-exp version
            this_alpha = 0.001 + (all_com(zzz,1) ./ (1+exp((90-v_ind)/(90./6)))) .* (exp(-v_ind/(all_com(zzz,2))));
        else
            % single exp version
            this_alpha = 0.001 + (all_com(zzz,1)*exp(-v_ind/(all_com(zzz,2))));
        end

        alpha = [alpha this_alpha];

        mouse(zz).opt_r2(zzz) = max(S.opt_r2_tensor,[],"all");
        mouse(zz).com_params(zzz,:) = all_com(zzz,:);

        % Quantify something like the total alpha for learning
        mouse(zz).tot_alpha(zzz) = trapz( alpha_da(mouse(zz).com_params(zzz,3),mouse(zz).com_params(zzz,1),mouse(zz).com_params(zzz,2),1:500) );

        % Get some insight into the optimal alpha
        [mouse(zz).opt_rew_poss(zzz),ind_rew_max]   = max(S.opt_RColl_tensor,[],"all"); 
        [a2_inds,a5_inds,a1_inds]                   = ind2sub(size(S.opt_r2_tensor),ind_rew_max);
        mouse(zz).opt_com_maxRew(zzz,:)             = [S.a2_vec(a2_inds), S.a5_vec(a5_inds), S.a1_vec(a1_inds)];
        mouse(zz).opt_rew_act(zzz)                  = S.tot_rew;

        disp(['Mouse: ' num2str(zz) ' | Session ' num2str(zzz) ': Model max rew: ' num2str(mouse(zz).opt_rew_poss(zzz)) ' ; Best fit rew: ' num2str(S.tot_rew)])

    end

    figure(900);
    subplot(6,6,unos(zz)); 
    yyaxis right; plot(times,alpha,'color',[0.9886    0.8096    0.1454]);
    
    %----------------------------------------
    %----------------------------------------
    % ------------- ALPHA fitting to DA
    figure(1000); clf;

    for sess=unique(rw_p_inds_s)'
            fitfun = fittype( alpha_da );
    
            targety = movmean(da_resp_rew(sess).trap,31)./max(movmean(da_resp_rew(1).trap,31));
    
            % reasonable initial guesses
            a0 = [ 0.2 0.5 500 ];
    
            [f,gof] = fit([1:numel(targety)]',targety,fitfun,'StartPoint',a0,'Upper',[1 1 2*numel(targety)],'Lower',[0 0 10]);  

            da_resp_rew(sess).alpha_f   = f;
            da_resp_rew(sess).r_sq      = gof.rsquare;

            subplot(1,5,sess);
            plot(f,[1:numel(targety)]',targety); box off;
            ylim([0 1.5]); xlim([0 500]);
            title(num2str(gof.rsquare));

            yyaxis right;
            plot(1:500,alpha_da(0.01,mouse(zz).com_params(sess,1),mouse(zz).com_params(sess,2),1:500),'b');
            ylim([0 1.5.*max(mouse(zz).com_params(:,1))]);

            % Quantify something like the total DA delivered to learning
            % (equivalent to total alpha delivered)
            mouse(zz).tot_dopa(sess) = trapz( alpha_da(f.a1,f.a2,f.a5,1:500) );
    end

    figure(901);
    subplot(6,2,zz);
    scatter(mouse(zz).tot_alpha,mouse(zz).tot_dopa,50,unique(rw_p_inds_s)','filled'); colormap(sess_map);
    axis([0 150 50 450]); box off;

    %----------------------------------------
    %----------------------------------------
    % Run new sims with optimized alpha over entire dataset
    %----------------------------------------
    %----------------------------------------
    num_iter            = 1; 
    vismat              = zeros(6,size(visit_matrix,2),num_iter);
    rewmat              = zeros(6,size(visit_matrix,2),num_iter);
    trans_r2_iter       = zeros(1,num_iter);
    income_r2_iter      = zeros(1,num_iter);
    income_model        = zeros(1,numel(income),num_iter);
    p_reward            = zeros(6,size(visit_matrix,2),num_iter);
    Q_reward            = zeros(6,size(visit_matrix,2),num_iter);

    cost_per_port =                 ...
    [1	14	18	70	72.2	65.5;   ...
    14	1	22.8	56	65.5	42; ...
    18	22.8	1	72.2	70	56; ...
    70	56	72.2	1	18	22.8;   ...
    72.2	65.5	70	18	1	14; ...
    65.5	42	56	22.8	14	1];

    alpha_Q = [mean(alpha) mean(alpha)]; 
    beta    = 1;

    for iter = 1:num_iter
        [trans_r2_iter(iter), income_r2_iter(iter), vismat(:,:,iter), rewmat(:,:,iter), p_reward(:,:,iter), income_model(:,:,iter)] = HX_model_session_forAlphaConcat(alpha,visit_matrix,cost_per_port.^mean(all_com(:,3)),rew_sched,income);
        [~, ~, ~, ~, Q_reward(:,:,iter)] = HX_model_session_QConcat(alpha_Q,beta,visit_matrix,cost_per_port,rew_sched,income);
    end

    visits_for_LL = squeeze(mean(vismat,3));
        visits_for_LL_sm = zeros(size(visits_for_LL));
    rewards_for_LL = squeeze(mean(rewmat,3));
        rewards_for_LL_sm = zeros(size(rewards_for_LL));
    p_reward_all = squeeze(mean(p_reward,3));
    Q_reward_all = squeeze(mean(Q_reward,3));

    figure(8); clf;
    plot(1:numel(all_visits),[0 diff(all_sessid)],'k-'); hold on;
    plot(1:numel(all_visits),income,'color',sess_map(3,:),'linewidth',2);
    plot(1:numel(all_visits),mean(income_model,3),'color',sess_map(end,:),'linewidth',2);
    xlim([0 numel(all_visits)]); box off;

    figure(900);
    subplot(6,6,tres(zz));
    plot(times,[0 diff(all_sessid)],'k-'); hold on;
    plot(times,income,'color',sess_map(3,:),'linewidth',2);
    plot(times,mean(income_model,3),'color',sess_map(end,:),'linewidth',2);
    xlim([0 5.4e4]); box off;

    uv_inds_m             = find(sum(visits_for_LL,1)==1);
    [~,ports_m]           = max(visits_for_LL(:,uv_inds_m),[],1);

    rw_inds_m             = find(sum(rewards_for_LL,1)==1);
    vis_are_rew           = ismember(uv_inds_m,rw_inds_m);
    rtimes_m              = times(vis_are_rew);
    [~,rports_m]          = max(rewards_for_LL(:,rw_inds_m),[],1);

    figure(11); clf;
    for qqq = 1:6
        visits_for_LL(qqq,ceil(times(ports_m==qqq)))=1;
        visits_for_LL_sm(qqq,:) = movmean(visits_for_LL(qqq,:),10*60);

        figure(11); subplot(211);
        plot(visits_for_LL_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 size(visits_for_LL,2)]); ylim([0 0.15]);
        ylabel('Fraction of visits');

        rewards_for_LL(qqq,ceil(rtimes_m(rports_m==qqq)))=1;
        rewards_for_LL_sm(qqq,:) = movmean(rewards_for_LL(qqq,:),30*60);

        figure(11); subplot(212);
        semilogy(rewards_for_LL_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 size(rewards_for_LL,2)]); ylim([0 0.075]);
        ylabel('Fraction of visits rewarded');
    end
    figure(11); subplot(211); yyaxis right; plot(session_ids,'k'); box off;
    ylabel('Session');
    all_rewards = sum(rewards_for_LL_sm,1);
    figure(11); subplot(212); plot(all_rewards,'k-'); yyaxis right; plot(session_ids,'k'); box off;
    ylabel('Session');

    aqua_model.rewards_for_LL_sm    = rewards_for_LL_sm;
    aqua_model.rewards_for_LL       = rewards_for_LL;
    aqua_model.visits_for_LL_sm     = visits_for_LL_sm;
    aqua_model.visits_for_LL        = visits_for_LL;   
    vis_inds = find(sum(squeeze(vismat(:,:,1)),1));
    aqua_model.rewards_invis        = rewards_for_LL(:,vis_inds);
    aqua_model.visits_invis         = visits_for_LL(:,vis_inds);    
    
    choice_predictions = mean( reshape( (visits_for_LL_sm-visit_matrix_sm).^2 , 1, numel(visits_for_LL_sm)) );
    null_predictions = zeros(1,1000);
    for qqq=1:1000
        tmp_inds = randperm(6);
        if sum(tmp_inds==[1 2 3 4 5 6])<6 
            null_predictions(qqq) = mean( reshape( (visits_for_LL_sm-visit_matrix_sm(tmp_inds,:)).^2 , 1, numel(visits_for_LL_sm)) );
        else
            tmp_inds=tmp_inds(6:-1:1);
            null_predictions(qqq) = mean( reshape( (visits_for_LL_sm-visit_matrix_sm(tmp_inds,:)).^2 , 1, numel(visits_for_LL_sm)) );
        end
    end

    sh_nom_pvalue = numel(find(null_predictions<choice_predictions)) ./ numel(null_predictions);

    %----------------------------------------
    %----------------------------------------
    % export the GLM data structure
    %----------------------------------------
    %----------------------------------------

    p_rew_s12           = squeeze(mean(p_reward(:,session_ids<3,:),3));
    Q_rew_s12           = squeeze(mean(Q_reward(:,session_ids<3,:),3));
    rw_times_s12        = T.hexa_data.event_time_con(rw_p_inds_s12);
 
    if numel(find(T.hexa_data.photo.sess==2))>0
        da_sink_s2          = TNC_ExtTrigWins(photo_data(T.hexa_data.photo.sess==2),rw_p_inds_i(rw_p_inds_s==2),photo_event_win);
    else
        da_sink_s2.inds = [];
        da_sink_s2.wins = [];
    end
    
    GLM_export(zz).anim_id      = all_sess_files(zz).name(1:brk-1);
    GLM_export(zz).all_com      = all_com;

    GLM_export(zz).alpha        = alpha(rw_p_logic_s12==1);
    GLM_export(zz).port_id      = T.hexa_data.port_n(rw_p_inds_s12);

    GLM_export(zz).aqua_model   = aqua_model;
    GLM_export(zz).loss_target  = loss_target;

    for qqq=1:numel(rw_times_s12)
        GLM_export(zz).AQUA_rpe(qqq) = 1-p_rew_s12(GLM_export(zz).port_id(qqq),round(rw_times_s12(qqq)));
        GLM_export(zz).Q_rpe(qqq)    = 1-Q_rew_s12(GLM_export(zz).port_id(qqq),round(rw_times_s12(qqq)));
    end

end

% save GLM_export_DMS_v7 GLM_export -v7
% save GLM_export_NAc_v7 GLM_export -v7

% save MOUSE_export_DMS_v7 mouse -v7
% save MOUSE_export_NAc_v7 mouse -v7

% save MOUSE_export_ALL_v7 mouse -v7