%% HX_final_script

%% Nice simple extraction of data into a convenient strcture concatenated over all sessions

%----------------------------------------
%----------------------------------------
% Let the data control which animal is being examined and then look over sessions for fit params within that    
%----------------------------------------
%----------------------------------------
all_sess_files = dir('*NAc*dat.mat');
% all_sess_files = dir('*DMS*dat.mat');
% all_sess_files = dir('*dat.mat');

figure(900); clf; figure(101); clf; unos = 1:3:34; dos = 2:3:35; tres = 3:3:36;
clear GLM_export mouse;

loss_target = 'combined'

for zz=1:numel(all_sess_files)

    %----------------------------------------
    %----------------------------------------
    % load the corresponding hexa_data_an     
    %----------------------------------------
    %----------------------------------------
    T = load(all_sess_files(zz).name);

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

    session_bounds      = find([0 diff(T.hexa_data.session_n')]==1);
    if numel(session_bounds)==0
        sess_start_times    = [];
    else
        sess_start_times    = round(T.hexa_data.event_time_con(session_bounds)');
    end

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
    % Get reward times with valid photometry
    %----------------------------------------
    %----------------------------------------
    rw_p_inds           = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==1 & ~isnan(T.hexa_data.photo_i));
    rw_p_logic          = ismember(rw_inds,rw_p_inds);
    
    ur_p_inds           = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==0 & ~isnan(T.hexa_data.photo_i));
    ur_p_logic          = ismember(uv_inds,ur_p_inds);

    rw_p_inds_i         = T.hexa_data.photo_i(rw_p_inds);
    rw_p_inds_s         = T.hexa_data.session_n(rw_p_inds);
    rw_p_inds_p         = T.hexa_data.port_n(rw_p_inds);

    ur_p_inds_i         = T.hexa_data.photo_i(ur_p_inds);
    ur_p_inds_s         = T.hexa_data.session_n(ur_p_inds);
    ur_p_inds_p         = T.hexa_data.port_n(ur_p_inds);

    %----------------------------------------
    %----------------------------------------
    % Pull photometry data
    %----------------------------------------
    %----------------------------------------
    photo_event_win     = [240 400];
    photo_data          = T.hexa_data.photo.dFF;
    highpass_photo      = 1;

    if highpass_photo
        d2 = designfilt("highpassiir",FilterOrder=12, ...
            HalfPowerFrequency=0.0001,DesignMethod="butter");
        yH = filtfilt(d2,photo_data);
        photo_data = yH;
    end

    figure(701); clf;
    for sess=unique(rw_p_inds_s)'
        figure(701);
        da_resp_rew(sess).da_sink= TNC_ExtTrigWins(photo_data(T.hexa_data.photo.sess==sess),rw_p_inds_i(rw_p_inds_s==sess),photo_event_win);
        da_resp_rew(sess).trap   = mean(da_resp_rew(sess).da_sink.wins(:,photo_event_win(1):photo_event_win(1)+120),2);
        % da_resp_rew(sess).trap   = mean(da_resp_rew(sess).da_sink.wins(:,photo_event_win(1):photo_event_win(1)+120),2)-mean(da_resp_rew(sess).da_sink.wins(:,photo_event_win(1)-65:photo_event_win(1)),2);
        da_resp_rew(sess).port   = rw_p_inds_p(rw_p_inds_s==sess);
        da_resp_rew(sess).visi   = rw_p_inds(rw_p_inds_s==sess);
        da_resp_rew(sess).vist   = T.hexa_data.event_time_con(da_resp_rew(sess).visi);
            subplot(1,5,1); plot(da_resp_rew(sess).da_sink.avg,'color',sess_map(sess,:), 'LineWidth',2); hold on; box off; ylabel('dF/F'); xlabel('Time from rew');
            subplot(1,5,2:5); scatter(da_resp_rew(sess).visi,da_resp_rew(sess).trap,25,rw_p_inds_p(rw_p_inds_s==sess),'filled','MarkerFaceAlpha',0.5); hold on;
            subplot(1,5,2:5); plot(da_resp_rew(sess).visi,movmean(da_resp_rew(sess).trap,31),'color',[sess_map(sess,:) 1],'LineWidth',2); hold on;
    
        da_resp_ure(sess).da_sink= TNC_ExtTrigWins(photo_data(T.hexa_data.photo.sess==sess),ur_p_inds_i(ur_p_inds_s==sess),photo_event_win);
        da_resp_ure(sess).trap   = mean(da_resp_ure(sess).da_sink.wins(:,photo_event_win(1):photo_event_win(1)+240),2)-mean(da_resp_ure(sess).da_sink.wins(:,1:30),2);
        da_resp_ure(sess).port   = ur_p_inds_p(ur_p_inds_s==sess);
        da_resp_ure(sess).visi   = ur_p_inds(ur_p_inds_s==sess);
        da_resp_ure(sess).vist   = T.hexa_data.event_time_con(da_resp_ure(sess).visi);
            subplot(1,5,1); plot(da_resp_ure(sess).da_sink.avg,'-', 'color',sess_map(sess,:)/2, 'LineWidth',2); hold on;
            subplot(1,5,2:5); scatter(da_resp_ure(sess).visi,da_resp_ure(sess).trap,10,ur_p_inds_p(ur_p_inds_s==sess),'filled','MarkerFaceAlpha',0.5); colormap(sess_map); hold on;
            subplot(1,5,2:5); plot(da_resp_ure(sess).visi,movmean(da_resp_ure(sess).trap,31),'color',[sess_map(sess,:) 1]/2,'LineWidth',2); hold on;
            ylabel('DA Response (Rew; Light) (UnRew; Dark)');
            xlabel('Session time (sec)');
            drawnow;

        figure(900);
            subplot(6,6,unos(zz));
            plot(da_resp_rew(sess).vist,movmean(da_resp_rew(sess).trap,31),'color',[sess_map(sess,:) 1],'LineWidth',2); hold on;
            plot(da_resp_ure(sess).vist,movmean(da_resp_ure(sess).trap,31),'color',[sess_map(sess,:) 1]/2,'LineWidth',2); hold on;
            xlabel('Session time (sec)'); xlim([0 180*5*60]);  box off; brk = strfind(all_sess_files(zz).name,'_'); title(all_sess_files(zz).name(1:brk(1)-1));
            if zz==1
                ylabel('DA Response (Rew; Light) (UnRew; Dark)');
            end
    end

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
    % Look at session 2->3 transition if it exists
    %----------------------------------------
    %----------------------------------------
    transition = find(session_ids==2 & [diff(session_ids) 0]==1);
    clear P Pm Q Qm;
    if numel(transition)>0 & ~strcmp('ML15',all_sess_files(zz).name(1:brk(1)-1))
        twin        = [ transition-1 10e3 ];
        KL_div_23   = zeros(6,numel([-twin(1):twin(2)]));
        KL_div_23m  = zeros(6,numel([-twin(1):twin(2)]));
        figure(101); 
        % clf; subplot(121);
        subplot(6,6,unos(zz));

        plot([0 0 0 0.08],'k-','linewidth',2); hold on;
        ylabel('P(choice)'); xlabel('Time from switch (s)'); title('Model');

        % subplot(122);
        subplot(6,6,dos(zz));

        plot([0 0 0 0.08],'k-','linewidth',2); hold on;
        ylabel('P(choice)'); xlabel('Time from switch (s)'); title('Data');

        for qqq = 1:6
            % subplot(121);
            subplot(6,6,unos(zz));
            plot(-twin(1):twin(2),visits_for_LL_sm(qqq,transition-twin(1):transition+twin(2)),'color',sess_map(qqq,:),'linewidth',1); hold on;
            box off;

            Pm = visits_for_LL_sm(qqq,transition-twin(1):transition+twin(2))+0.0001;
            Qm = mean(visits_for_LL_sm(qqq,transition-4.5e3:transition));
            KL_div_23m(qqq,:) = Pm.*log(Pm/Qm);            
            
            % subplot(122);
            subplot(6,6,dos(zz));
            plot(-twin(1):twin(2),visit_matrix_sm(qqq,transition-twin(1):transition+twin(2)),'color',sess_map(qqq,:),'linewidth',1); hold on;
            box off;

            P = visit_matrix_sm(qqq,transition-twin(1):transition+twin(2))+0.0001;
            Q = mean(visit_matrix_sm(qqq,transition-4.5e3:transition));
            KL_div_23(qqq,:) = P.*log(P/Q);            
        end

        subplot(6,6,unos(zz));
        yyaxis right;
        plot(-twin(1):twin(2),sum(KL_div_23m,1),'color',[1 0 0 0.5],'LineWidth',3);
        
        subplot(6,6,dos(zz));
        yyaxis right;
        plot(-twin(1):twin(2),sum(KL_div_23,1),'color',[1 0 0 0.5],'LineWidth',3);
        

        subplot(6,6,tres(zz));
        plot(da_resp_rew(1).vist,movmean(da_resp_rew(1).trap,31),'color',[sess_map(1,:) 1],'LineWidth',3); hold on;
        plot(da_resp_rew(2).vist,movmean(da_resp_rew(2).trap,31),'color',[sess_map(2,:) 1],'LineWidth',3); hold on;
        plot(da_resp_rew(3).vist,movmean(da_resp_rew(3).trap,31),'color',[sess_map(3,:) 1],'LineWidth',3); hold on;
        ylim([-0.5 3]); title(all_sess_files(zz).name(1:brk(1)-1)); box off;

    end

    %----------------------------------------
    %----------------------------------------
    % export the GLM data structure
    %----------------------------------------
    %----------------------------------------

    av_inds_12          = find(T.hexa_data.unique_vis==1 & T.hexa_data.session_n<3);
    rw_p_inds_s12       = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==1 & T.hexa_data.session_n<3 & ~isnan(T.hexa_data.photo_i));
    rw_p_logic_s12      = ismember(av_inds_12,rw_p_inds_s12);
    p_rew_s12           = squeeze(mean(p_reward(:,session_ids<3,:),3));
    Q_rew_s12           = squeeze(mean(Q_reward(:,session_ids<3,:),3));
    
    rw_times_s12        = T.hexa_data.event_time_con(rw_p_inds_s12);
    da_sink_s1          = TNC_ExtTrigWins(photo_data(T.hexa_data.photo.sess==1),rw_p_inds_i(rw_p_inds_s==1),photo_event_win);
    if numel(find(T.hexa_data.photo.sess==2))>0
        da_sink_s2          = TNC_ExtTrigWins(photo_data(T.hexa_data.photo.sess==2),rw_p_inds_i(rw_p_inds_s==2),photo_event_win);
    else
        da_sink_s2.inds = [];
        da_sink_s2.wins = [];
    end
    
    GLM_export(zz).uv_inds      = [da_sink_s1.inds ; da_sink_s2.inds];
    GLM_export(zz).DA_waveforms = [da_sink_s1.wins ; da_sink_s2.wins];
    GLM_export(zz).DA_samples   = da_sink_s1.range;
    GLM_export(zz).anim_id      = all_sess_files(zz).name(1:brk-1);
    GLM_export(zz).all_com      = all_com;

    GLM_export(zz).da_resp_rew  = da_resp_rew;
    GLM_export(zz).da_resp_ure  = da_resp_ure;
    
    GLM_export(zz).alpha        = alpha(rw_p_logic_s12==1);
    GLM_export(zz).alpha_da     = alpha_da;
    GLM_export(zz).port_id      = T.hexa_data.port_n(rw_p_inds_s12);

    GLM_export(zz).aqua_model   = aqua_model;
    GLM_export(zz).loss_target  = loss_target;

    GLM_export(zz).KL_div_23m   = KL_div_23m;
    GLM_export(zz).KL_div_23    = KL_div_23;

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

%% Da response magnitude & Inferred alpha & P(choice|rank) smoothed and averaged across mice (focus esp. session 2 and 3)

figure(950); clf;
figure(949); clf;
total_da    = zeros(numel(GLM_export),5);
da_wf       = NaN.*ones(numel(GLM_export),numel(GLM_export(zz).da_resp_rew(1).da_sink.range),5);
da_int      = NaN.*ones(numel(GLM_export),5);
anim_map = TNC_CreateRBColormap(numel(GLM_export),'shadowplay');
all_data.da_trap = zeros(numel(GLM_export),5.5e4);

for zz=1:numel(GLM_export)

    all_data.da_trap(zz,1:round(GLM_export(zz).da_resp_rew(1).vist(2))) = NaN;

    for sess = 1:numel(GLM_export(zz).da_resp_rew)
        if strcmp('ML15',GLM_export(zz).anim_id) & sess==3
            all_data.da_trap(zz,round(GLM_export(zz).da_resp_rew(sess).vist(1)):round(GLM_export(zz).da_resp_rew(sess).vist(end))) = NaN;            
            total_da(zz,sess) = NaN;
        else
            all_data.da_trap(zz,round(GLM_export(zz).da_resp_rew(sess).vist)) = movmean(GLM_export(zz).da_resp_rew(sess).trap,31);
            if numel(GLM_export(zz).da_resp_rew(sess).vist)>0
                total_da(zz,sess) = mean(GLM_export(zz).da_resp_rew(sess).trap);
            end

            da_wf(zz,:,sess) = GLM_export(zz).da_resp_rew(sess).da_sink.avg;
            da_int(zz,sess) = trapz(GLM_export(zz).da_resp_rew(sess).da_sink.avg(1,250:510));
            figure(949); 
            subplot(1,5,sess);
            plot(GLM_export(zz).da_resp_rew(sess).da_sink.range,GLM_export(zz).da_resp_rew(sess).da_sink.avg); hold on;
            axis([GLM_export(zz).da_resp_rew(sess).da_sink.range(1) GLM_export(zz).da_resp_rew(sess).da_sink.range(end) -1 3]);
        end

    end

    total_da(total_da==0) = NaN;

    last_t = max( round(GLM_export(zz).da_resp_rew(sess).vist) );

    % fill in missing session values with NaNs
    all_data.da_trap(zz,last_t:end) = NaN;

    % interpolate between known values
    known_x = find(all_data.da_trap(zz,1:last_t)~=0);
    known_y = all_data.da_trap(zz,known_x);

    unknown_x = find(all_data.da_trap(zz,1:last_t)==0);
    unknown_y = interp1(known_x,known_y,unknown_x);

    all_data.da_trap(zz,unknown_x) = unknown_y;

    figure(950); 
    plot(1:last_t,all_data.da_trap(zz,1:last_t),'color',[anim_map(zz,:) 0.33],'LineWidth',1); hold on;
    % plot(1:last_t,all_data.da_trap(zz,1:last_t),'color',[0 0 0 0.25],'LineWidth',1); hold on;
    % plot(da_resp_ure(sess).vist,movmean(da_resp_ure(sess).trap,31),'color',[sess_map(sess,:) 1]/2,'LineWidth',2); hold on;
    xlabel('Session time (sec)'); xlim([0 180*5*60]);  box off;

end

figure(950);
plot(mean(all_data.da_trap,1,'omitnan'),'color',[0 0 0],'LineWidth',5);

figure(948); clf;

for sess=1:5
    figure(949); 
    subplot(1,5,sess);
    plot(GLM_export(zz).da_resp_rew(sess).da_sink.range,mean(da_wf(:,:,sess),1,'omitnan'),'color',[0 0 0],'LineWidth',5); box off; hold on;
    ylim([-0.5 2.5]);

    if sess==2 | sess==3
        da_transients(sess,:) = squeeze(mean(da_wf(:,:,sess),1,'omitnan'));
        figure(948);
        plot(GLM_export(zz).da_resp_rew(sess).da_sink.range,mean(da_wf(:,:,sess),1,'omitnan'),'color',sess_map(sess,:),'LineWidth',5); hold on;
    end
end


valids = find(~isnan(total_da(:,3)))
[h,p] = ttest2(total_da(valids,2),total_da(valids,3))

figure(951); clf;
boxchart(total_da);

figure(952); clf;
boxchart(da_int)
[h,p] = ttest2(da_int(valids,2),da_int(valids,3))

%% Look at KL_div around transition and compare to dopamine directly

figure(970); clf;
figure(971); clf;
cnt = 1;

for zz=1:numel(GLM_export)

    if numel(GLM_export(zz).da_resp_rew)>2 & ~strcmp('ML15',GLM_export(zz).anim_id)
        
        transition      = max(GLM_export(zz).da_resp_rew(2).vist);
        transition2     = max(GLM_export(zz).da_resp_rew(1).vist);
        time_around     = 3000;

        figure(970);
        subplot(8,1,cnt);
        plot(GLM_export(zz).da_resp_rew(1).vist,movmean(GLM_export(zz).da_resp_rew(1).trap,31),'color',[sess_map(1,:) 1],'LineWidth',2); hold on;
        plot(GLM_export(zz).da_resp_rew(2).vist,movmean(GLM_export(zz).da_resp_rew(2).trap,31),'color',[sess_map(2,:) 1],'LineWidth',2); hold on;
        plot(GLM_export(zz).da_resp_rew(3).vist,movmean(GLM_export(zz).da_resp_rew(3).trap,31),'color',[sess_map(3,:) 1],'LineWidth',2); hold on;
        ylim([0 3.3]);

        yyaxis right;
        plot(sum(GLM_export(zz).KL_div_23,1),'color',[1 0 0 0.5],'LineWidth',2); 
        ylim([-0.1 0.4]); box off;

        xlim([transition-time_around transition+time_around]); title(transition);
        drawnow;

        stable_1_ind    = find(GLM_export(zz).da_resp_rew(1).vist>transition2-time_around,1,'first');
        init_2_ind      = find(GLM_export(zz).da_resp_rew(2).vist-GLM_export(zz).da_resp_rew(2).vist(1)>time_around,1,'first');
        
        stable_2_ind    = find(GLM_export(zz).da_resp_rew(2).vist>transition-time_around,1,'first');
        init_3_ind      = find(GLM_export(zz).da_resp_rew(3).vist-GLM_export(zz).da_resp_rew(3).vist(1)>time_around,1,'first');

        da_end_2        = mean(GLM_export(zz).da_resp_rew(2).trap(stable_2_ind:end));
        da_start_3      = mean(GLM_export(zz).da_resp_rew(3).trap(1:init_3_ind));
        da_end_1        = mean(GLM_export(zz).da_resp_rew(1).trap(stable_1_ind:end));
        da_start_2      = mean(GLM_export(zz).da_resp_rew(2).trap(1:init_2_ind));

        kl_end_2        = mean(sum(GLM_export(zz).KL_div_23(:,transition-time_around:transition),1));
        kl_start_3      = mean(sum(GLM_export(zz).KL_div_23(:,transition:transition+time_around),1));
        kl_end_1        = mean(sum(GLM_export(zz).KL_div_23(:,transition2-time_around:transition2),1));
        kl_start_2      = mean(sum(GLM_export(zz).KL_div_23(:,transition2:transition2+time_around),1));

        what_weve_got(cnt,1) = da_end_1;
        what_weve_got(cnt,2) = da_start_2;
        what_weve_got(cnt,3) = da_end_2;
        what_weve_got(cnt,4) = da_start_3;

        what_weve_got(cnt,5) = kl_end_1;
        what_weve_got(cnt,6) = kl_start_2;
        what_weve_got(cnt,7) = kl_end_2;
        what_weve_got(cnt,8) = kl_start_3;
        
        figure(971);
        % subplot(8,1,cnt);
        plot([da_start_2 da_end_2 da_start_3],[ kl_start_2 kl_end_2 kl_start_3],'color',[0.5 0.5 0.5 0.1],'LineWidth',3); hold on;
        scatter([da_start_2 da_end_2 da_start_3],[ kl_start_2 kl_end_2 kl_start_3],100,sess_map(2:4,:),'filled'); hold on;
        if numel(strfind(GLM_export(zz).anim_id,'ML'))>0
            scatter([ da_start_2 da_end_2 da_start_3],[ kl_start_2 kl_end_2 kl_start_3],100,'k'); hold on;
        end

        

        cnt = cnt +1;
    end
end

% [R,P]=corrcoef([what_weve_got(:,2)],[what_weve_got(:,4)]);
[R,P]=corrcoef([what_weve_got(:,1) ; what_weve_got(:,2) ; what_weve_got(:,3) ; what_weve_got(:,4)],[what_weve_got(:,5) ; what_weve_got(:,6) ; what_weve_got(:,7) ; what_weve_got(:,8)])
xlabel('DA response'); ylabel('KL Divergence from session 2 stable (s.d.)');
title(['rho= ' num2str(R(1,2)) ' ; p=' num2str(P(1,2))])

da_and_kl_dat_nac = what_weve_got;
save ForLauraFinalFig da_and_kl_dat_nac -v7

% kruskalwallis(what_weve_got(:,[1 2]))

%% Run a simulated sessions data with varying alpha to see the "ground truth" of how KL divergence should be related to alpha in principle
check_size = [0.25 0.5 0.75 1 1.5 2];
    figure(1051); clf;
    figure(1021); clf;

    clear store_visits;

        time_around     = 900;

for rep = 1:5
for cc=1:6

    alpha_check = alpha.*check_size(cc);

    for iter = 1:num_iter
        [trans_r2_iter(iter), income_r2_iter(iter), vismat(:,:,iter), rewmat(:,:,iter), p_reward(:,:,iter), income_model(:,:,iter)] = HX_model_session_forAlphaConcat(alpha_check + (rand(1)./10),visit_matrix,cost_per_port.^mean(all_com(:,3)),rew_sched,income);
    end

    visits_for_LL = squeeze(mean(vismat,3));
        visits_for_LL_sm = zeros(size(visits_for_LL));
    rewards_for_LL = squeeze(mean(rewmat,3));
        rewards_for_LL_sm = zeros(size(rewards_for_LL));

    figure(1008); clf;
    plot(1:numel(all_visits),[0 diff(all_sessid)],'k-'); hold on;
    plot(1:numel(all_visits),income,'color',sess_map(3,:),'linewidth',2);
    plot(1:numel(all_visits),mean(income_model,3),'color',sess_map(end,:),'linewidth',2);
    xlim([0 numel(all_visits)]); box off;

    uv_inds_m             = find(sum(visits_for_LL,1)==1);
    [~,ports_m]           = max(visits_for_LL(:,uv_inds_m),[],1);

    rw_inds_m             = find(sum(rewards_for_LL,1)==1);
    vis_are_rew           = ismember(uv_inds_m,rw_inds_m);
    rtimes_m              = times(vis_are_rew);
    [~,rports_m]          = max(rewards_for_LL(:,rw_inds_m),[],1);

    figure(1011); clf;
    for qqq = 1:6
        visits_for_LL(qqq,ceil(times(ports_m==qqq)))=1;
        visits_for_LL_sm(qqq,:) = movmean(visits_for_LL(qqq,:),10*60);

        figure(1011); subplot(211);
        plot(visits_for_LL_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 size(visits_for_LL,2)]); ylim([0 0.15]);
        ylabel('Fraction of visits');

        rewards_for_LL(qqq,ceil(rtimes_m(rports_m==qqq)))=1;
        rewards_for_LL_sm(qqq,:) = movmean(rewards_for_LL(qqq,:),30*60);

        figure(1011); subplot(212);
        semilogy(rewards_for_LL_sm(qqq,:),'color',sess_map(qqq,:),'linewidth',2); hold on;
        xlim([0 size(rewards_for_LL,2)]); ylim([0 0.075]);
        ylabel('Fraction of visits rewarded');
    end
    all_rewards = sum(rewards_for_LL_sm,1);

    store_visits(cc).vis_sm(:,:,rep) = visits_for_LL_sm;

    %----------------------------------------
    %----------------------------------------
    % Look at session 2->3 transition if it exists
    %----------------------------------------
    %----------------------------------------
    transition = find(session_ids==2 & [diff(session_ids) 0]==1);
    transition2 = find(session_ids==1 & [diff(session_ids) 0]==1);
    
    clear Pm Qm;
    if numel(transition)>0 & ~strcmp('ML15',all_sess_files(zz).name(1:brk(1)-1))
        twin        = [ transition-1 10e3 ];
        KL_div_23m  = zeros(6,numel([-twin(1):twin(2)]));

        for qqq = 1:6
            Pm = visits_for_LL_sm(qqq,transition-twin(1):transition+twin(2))+0.0001;
            Qm = mean(visits_for_LL_sm(qqq,transition-4.5e3:transition));
            KL_div_23m(qqq,:) = Pm.*log(Pm/Qm);                        
        end
    end
    
    figure(1021);
    plot(sum(KL_div_23m,1),'linewidth',2,'color',[sess_map(cc,:) 0.8]); hold on;
    % xlim([transition-2000 transition+2000]); 
    box off;
    ylabel('KL Divergence'); xlabel('Session time');

    klMC_end_2(cc,rep)        = mean(sum(KL_div_23m(:,transition-time_around:transition),1));
    klMC_start_3(cc,rep)      = mean(sum(KL_div_23m(:,transition:transition+time_around),1));

    klMC_end_1(cc,rep)        = mean(sum(KL_div_23m(:,transition2-time_around:transition2),1));
    klMC_start_2(cc,rep)      = mean(sum(KL_div_23m(:,transition2:transition2+time_around),1));
    
    vis_ts          = find(sum(visit_matrix,1)==1);
    vis_trans       = find(vis_ts>transition,1)-1;
    vis_trans2      = find(vis_ts>transition2,1)-1;
    vis_taround     = 150;
    
    alph_end_1(cc,rep)        = mean(alpha_check(vis_trans2-vis_taround:vis_trans2));
    alph_start_2(cc,rep)      = mean(alpha_check(vis_trans2:vis_trans2+vis_taround));

    alph_end_2(cc,rep)        = mean(alpha_check(vis_trans-vis_taround:vis_trans));
    alph_start_3(cc,rep)      = mean(alpha_check(vis_trans:vis_trans+vis_taround));
    
    figure(1051);
    plot([alph_start_2(cc,rep) alph_end_2(cc,rep) alph_start_3(cc,rep)],[klMC_start_2(cc,rep) klMC_end_2(cc,rep) klMC_start_3(cc,rep)],'color',[0.5 0.5 0.5 0.1],'LineWidth',1); hold on;

end

    scatter(alph_end_2(:,rep),klMC_end_2(:,rep),50,sess_map(3,:),'filled'); hold on;
    scatter(alph_start_3(:,rep),klMC_start_3(:,rep),50,sess_map(4,:),'filled');
    % scatter(alph_end_1(:,rep),klMC_end_1(:,rep),50,sess_map(1,:),'filled');
    scatter(alph_start_2(:,rep),klMC_start_2(:,rep),50,sess_map(2,:),'filled');
   
    ylabel('KL Divergence (std)'); xlabel('Nominal \alpha'); box off;

end

alpha_and_kl_6pg5_aqua.alph_start_2 = alph_start_2;
alpha_and_kl_6pg5_aqua.alph_end_2 = alph_end_2;
alpha_and_kl_6pg5_aqua.alph_start_3 = alph_start_3;

alpha_and_kl_6pg5_aqua.klMC_start_2 = klMC_start_2;
alpha_and_kl_6pg5_aqua.klMC_end_2 = klMC_end_2;
alpha_and_kl_6pg5_aqua.klMC_start_3 = klMC_start_3;

save ForLauraFinalFig2 alpha_and_kl_6pg5_aqua -v7


%% --
% Examine single port transitions over complete sessions for final figure
% use per_sess_rank & store_visits
figure(670); clf;
for cc=1:6
    subplot(1,6,cc);
    shadedErrorBar(1:size(store_visits(cc).vis_sm,2),mean(store_visits(cc).vis_sm(cc,:,:),3),std(store_visits(cc).vis_sm(cc,:,:),[],3));   
    hold on;
    plot(1:size(visit_matrix_sm,2),visit_matrix_sm(cc,:),'r');
    ylim([0 0.1]); box off;
end


%% Take saved files and calculate summary stats for session 2 (r2 across a range of models)

%--------------
%--------------
% SESSION TO ANALYZE
sess=1  ;
%--------------
%--------------

a=1;
b=1;
c=0;
frac = 0.99;
clear opt
sym = TNC_CreateRBColormap(8,'bb-sym');

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

    loss = a*S.opt_r2_tensor-b*S.opt_inc_tensor;
    
    [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));
    sym = TNC_CreateRBColormap(numel(top_xperc_inds),'cpb');
        
    [a2_inds,a5_inds,a1_inds] = ind2sub(size(S.opt_r2_tensor),top_xperc_inds);
    
    [com_in_param_space] = centerOfMass3D(S.a2_vec(a2_inds), S.a5_vec(a5_inds), S.a1_vec(a1_inds), loss(top_xperc_inds)');

    figure(799+sess); clf;
    subplot(121);
    hold on;
    scatter3(com_in_param_space(3),com_in_param_space(2),com_in_param_space(1),100,mean(S.opt_r2_tensor(top_xperc_inds)'),'filled'); grid on;
    axis([min(S.a1_vec) max(S.a1_vec) min(S.a5_vec) max(S.a5_vec) min(S.a2_vec) max(S.a2_vec) ]); colormap(sym); colorbar;
    xlabel('a1'); ylabel('a5'); zlabel('a2'); view([-40 15]);

    subplot(122);
    x = 1:1000;
    alpha_vis = com_in_param_space(3) + (com_in_param_space(1)*exp(-x/com_in_param_space(2)));    
    where_r2 = 2*(max(S.opt_r2_tensor(top_xperc_inds)')-0.3);
    if where_r2>1
        where_r2=1;
    elseif where_r2<0
        where_r2=0;
    end
    plot(1:1000,alpha_vis,'color',[where_r2,0.5.*(1-where_r2),(1-where_r2)]); hold on;
    axis([0 1000 0 0.33]);
    ylabel('alpha(rew)'); xlabel('Num rewards'); box off;
    title([all_sess_files(zz).name(1:breaks(1)-1) ' ' all_sess_files(zz).name(breaks(2)-2:breaks(2)-1)]);
    drawnow; pause(0.1);

    all_r2      = [all_r2 ; max(S.opt_r2_tensor,[],"all")];
    all_coms    = [all_coms ; com_in_param_space];
    all_isos    = [all_isos ; numel(top_xperc_inds)];

    %------------------------------------------------
    % NEED TO ADD IN COMPARISON OF INCOME MSE FOR OPTIMAL COMPARED TO BEST/COMPARABLE FIXED ALPHA
opt.inc_diff_dyn(zz,1)  = min(S.opt_inc_tensor,[],"all");
opt.loss_diff_dyn(zz,1)  = max(loss,[],"all");

opt.inc_diff_flat(zz,1) = min(S.opt_inc_tensor(1,:,:),[],"all");
opt.loss_diff_flat(zz,1) = max(loss(1,:,:),[],"all");

[~,top_ind] = max(S.opt_inc_tensor,[],"all");        
[a2_ind,a5_ind,a1_ind] = ind2sub(size(S.opt_inc_tensor),top_ind);
opt.inc_diff_flat2(zz,1) = S.opt_inc_tensor(1, a5_ind, a1_ind);
opt.loss_diff_flat2(zz,1) = loss(1, a5_ind, a1_ind);

    % load the corresponding hexa_data_an struct
    Z = load([all_sess_files(zz).name(1:end-7) 'an.mat']);

    % Get data/optimal discovered fit observations for comparison
    opt.r2(zz,1) = max(S.opt_r2_tensor,[],"all");
    opt.inc(zz,1) = min(S.opt_inc_tensor,[],"all");
    opt.labels{1} = 'AQUA opt grid';
    opt.rew(zz,1) = median(S.opt_RColl_tensor(top_xperc_inds));
    opt.rew_act(zz,1) = numel(find(Z.hexa_data_an.rewards==1));
    opt.rew_aqua_opt(zz,1) = max(S.opt_RColl_tensor,[],"all");


        %----------------------------------------
        %-------------- FITTING ALPHA TO DA
        %----------------------------------------
        alpha = @(a1,a2,a5,x) a1 + (a2*exp(-x/a5));
        fitfun = fittype( alpha );        
        % data = movmean(Z.hexa_data_an.da_resp_all.r,25);
        data = Z.hexa_data_an.da_resp_all.r;
        targety = data./mean(data(1:75));
        targetx = find(Z.hexa_data_an.da_resp_data.r_vis_id==1);
        if numel(targetx)<600
            end_pt = numel(targetx);
        else
            end_pt = 750;
        end
        a0 = [ 0 1 100 ];        
        [f,gof] = fit(targetx(1:end_pt),targety(1:end_pt),fitfun,'StartPoint',a0,'Upper',[0.4 0.6 750],'Lower',[0.01 0.01 5])
        %----------------------------------------
        %-------------- FITTING ALPHA TO DA
        %----------------------------------------


% ------------- 
% ------------- RUN RANGE OF MODEL SIMULATIONS

        num_iter = 20;

        all_visits = find(sum(Z.hexa_data_an.visits,1)==1);
        rew_logic = sum(Z.hexa_data_an.rewards,1);
        all_rewards = rew_logic(all_visits);

        vismat              = zeros(6,numel(all_visits),num_iter);
        rewmat              = zeros(6,numel(all_visits),num_iter);
        trans_r2_iter       = zeros(1,num_iter);
        income_r2_iter      = zeros(1,num_iter);
        p_reward            = zeros(6,size(Z.hexa_data_an.visits,2),num_iter);

        [~,close_a2] = min(abs(f.a2-S.a2_vec))
        [~,close_a5] = min(abs(f.a5-S.a5_vec))

        prior = S.prior;
        
        [~,a1_ind_DA] = max(S.opt_r2_tensor(close_a2,close_a5,:));
        best_a1 = S.a1_vec(a1_ind_DA);
        best_a2 = S.a2_vec(close_a2);
        best_a5 = S.a5_vec(close_a5);

        % run AQUA for optimal alpha fit
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter), p_reward(:,:,iter)] = HX_model_session_forAlphaOpt(com_in_param_space(3),com_in_param_space(1),com_in_param_space(1),0,com_in_param_space(2),'sig_exp',S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,2) = median(trans_r2_iter);
        opt.rank(zz,1) = numel(find(S.opt_r2_tensor<=median(trans_r2_iter)))./prod(size(S.opt_r2_tensor)); % percentile
        opt.labels{2} = 'AQUA opt com';
        opt.rew(zz,2) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,2) = median(income_r2_iter);

        % FOR LAURA I NEED TO SAVE RPE (1-P_rew) per reward response from dopamine
        % save a separate file of csv data that she can use for GLM building

        GLM_export(zz).da_traces            = Z.hexa_data_an.da_resp_data.wins(Z.hexa_data_an.da_resp_data.r_vis_id==1,:);

            sample_logic                    = sum(S.visit_matrix,1);
            tmp_rewards                     = sum(Z.hexa_data_an.rewards,1);
            rew_logic                       = tmp_rewards(sample_logic==1);
            v_ind                           = 1:sum(sample_logic);    
            alpha_func                      = alpha(com_in_param_space(3),com_in_param_space(1),com_in_param_space(2),v_ind);
        GLM_export(zz).alpha                = alpha_func(rew_logic==1)';
        GLM_export(zz).alpha                = GLM_export(zz).alpha(1:size(GLM_export(zz).da_traces,1));

            rpe                             = 1-mean(p_reward,3);
        GLM_export(zz).rpe                  = rpe(Z.hexa_data_an.rewards==1 & S.visit_matrix==1);
        GLM_export(zz).rpe                  = GLM_export(zz).rpe(1:size(GLM_export(zz).da_traces,1));
        
            breaks                          = strfind(all_sess_files(zz).name,'_');
        GLM_export(zz).name                 = all_sess_files(zz).name(1:breaks(2)-1);
        GLM_export(zz)

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
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(S.a1_vec(end).*f.a1,S.a2_vec(end).*f.a2,S.a2_vec(end).*f.a2,0,f.a5,'sig_exp',S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,3) = median(trans_r2_iter);
        opt.labels{3} = 'AQUA DA=alpha';
        opt.rew(zz,3) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,3) = mean(income_r2_iter);

        % run AQUA with flat alpha and no distance
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(com_in_param_space(3),0,0,0,0,'sig_exp',S.visit_matrix,ones(size(S.cost_per_port)),S.rew_sched,S.income,S.prior);
        end

        opt.r2(zz,4) = median(trans_r2_iter);
        opt.labels{4} = 'AQUA -dist -dyn';
        opt.rew(zz,4) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,4) = mean(income_r2_iter);
        
        % run AQUA without distance weighting
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_forAlphaOpt(com_in_param_space(3),com_in_param_space(1),com_in_param_space(1),0,com_in_param_space(2),'sig_exp',S.visit_matrix,ones(size(S.cost_per_port)),S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,5) = median(trans_r2_iter);
        opt.labels{5} = 'AQUA -dist';
        opt.rew(zz,5) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,5) = mean(income_r2_iter);

        % run Q learn static alpha optimum
        alpha = [com_in_param_space(3) com_in_param_space(3)]; 
        beta = 1; % Softmax term 
        % Grossman model uses ~ Beta = 4, alpha(1) = 0.5, alpha(2) = 0.1
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_Q(alpha,beta,S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,6) = median(trans_r2_iter);
        opt.labels{6} = 'Q-learn';
        opt.q.params.alpha = alpha; 
        opt.q.params.beta = beta; 
        opt.rew(zz,6) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,6) = mean(income_r2_iter);

        % run WSLS
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_WSLS(S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,7) = median(trans_r2_iter);
        opt.labels{7} = 'WSLS';
        opt.rew(zz,7) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,7) = mean(income_r2_iter);
        
        % run OIO
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_OIO(S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,8) = median(trans_r2_iter);
        opt.labels{8} = 'OIO';
        opt.rew(zz,8) = numel(find(rewmat==1))./num_iter;
        opt.inc(zz,8) = mean(income_r2_iter);

        
        % run RANDOM
        parfor iter = 1:num_iter
            [trans_r2_iter(1,iter),income_r2_iter(1,iter), vismat(:,:,iter),rewmat(:,:,iter)] = HX_model_session_RAND(S.visit_matrix,S.cost_per_port,S.rew_sched,S.income,S.prior);
        end
        opt.r2(zz,9) = median(trans_r2_iter);
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

            % plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',targety,'.','color', [0.8 0.8 0.8] );
            % hold on;
            % plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',movmean(targety,21), 'k-')
            % plot([1:numel(Z.hexa_data_an.da_resp_all.r)]',f([1:numel(Z.hexa_data_an.da_resp_all.r)]'),'r','linewidth',3);
            

            plot(targetx,targety,'.','color', [0.8 0.8 0.8] );
            hold on;
            plot(targetx,movmean(targety,21), 'k-')
            plot(targetx,f([1:numel(Z.hexa_data_an.da_resp_all.r)]'),'r','linewidth',2);
            
            % plot(f,x(round(end/20):end),y(round(end/20):end)); title(num2str(-1/f.d)); 
            ylim([-0.25 1.25]); title(['a1: ' num2str(f.a1) ' a2: ' num2str(f.a2) ' a5: ' num2str(f.a5)]);

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

clear session
session(sess).opt       = opt;
session(sess).all_coms  = all_coms;
session(sess).all_taus  = all_taus;
session(sess).all_r2    = all_r2;
session(sess).all_recloc= all_recloc;


save SessionDataCompareModels session -v7

% Look at boxchart for the new runs

figure((sess*100)+57); clf;
boxchart(opt.r2(all_recloc==1,[1:2 6:9]).^2);
xticklabels(opt.labels([1:2 6:9]));
ylabel('Transition matrix similarity (r^2)');
xlabel('Model type');
ylim([-0.1 1]);
title(sess);

[p,t,stats] = anova1(opt.r2(all_recloc==1,[1:2 6:9]));
[c,m,h,~] = multcompare(stats);

figure((sess*100)+58); clf;
boxchart([opt.rew_aqua_opt opt.rew(:,[2 6:9])]-opt.rew_act);
xticklabels(['AQUA optimium' opt.labels([2 6:9])]);
ylabel('\Delta Predicted Rewards Collected');
xlabel('Model type');
ylim([-200 100]);
title(sess);

[p_rew,t_rew,stats_rew] = anova1([opt.rew_aqua_opt opt.rew(:,[2 6:9])]-opt.rew_act);
[c_rew,m_rew,h_rew,~] = multcompare(stats_rew);

figure((sess*100)+59); clf;
subplot(121);
boxchart(opt.r2(all_recloc==1,2:3).^2);
xticklabels(opt.labels([2:3]));
ylabel('Transition matrix similarity (r^2)');
ylim([0 1]);
title(sess);

subplot(122);
boxchart(opt.rew(all_recloc==1,[2:3])-opt.rew_act(all_recloc==1));
xticklabels(opt.labels([2:3]));
ylabel('\Delta Predicted Rewards Collected');
ylim([-200 100]);


figure((sess*100)+60); clf;
boxchart(opt.r2(:,[2 5 4]).^2);
xticklabels(opt.labels([2 5 4]));
ylabel('Transition matrix similarity (r^2)');
xlabel('Model type');
ylim([-0.1 1]);
title(sess);

[p,t,stats_aq] = anova1(opt.r2(:,[2 5 4]).^2);
[c_aqua,m_aqua,h_aqua,~] = multcompare(stats_aq);

sess_cmap = TNC_CreateRBColormap(8,'hue7');
figure(161); 
if sess==1
    clf;
end
% subplot(121);
% plot([0.033 0.133], [0.033 0.133]); hold on;
% scatter(opt.inc_diff_dyn(:,1), opt.inc_diff_flat2(:,1),'filled');
% % xticklabels({'Dynamic','Constant'});
% % ylabel('Income RMSE');
% % xlabel('Alpha type');
% % ylim([-0.1 1]);
% title(sess);
% subplot(122);
plot([-0.25 1], [-0.25 1],'k--'); hold on;
% scatter([zeros(size(opt.loss_diff_dyn(:,1),1),1); ones(size(opt.loss_diff_dyn(:,1),1),1)],[opt.loss_diff_dyn(:,1); opt.loss_diff_flat2(:,1)],'k');
scatter(opt.loss_diff_dyn(:,1), opt.loss_diff_flat2(:,1),75,sess_cmap(sess,:),'filled');
axis([-0.25 1 -0.25 1]);
xlabel('Objective (dynamic)'); ylabel('Objective (static)');
% title(sess);

session(sess).p_dyn_sign = signrank(opt.inc_diff_flat2,opt.inc_diff_dyn);
session(sess).p_dyn_ranks = ranksum(opt.inc_diff_flat2,opt.inc_diff_dyn);
session(sess).opt = opt;
session(sess).c_r2 = c;
session(sess).c_rew = c_rew;
session(sess).c_aqua = c_aqua;