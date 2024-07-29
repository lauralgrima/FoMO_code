%% HX_final_script

%% Nice simple extraction of data into a convenient strcture concatenated over all sessions
% think about modeling the whole concatenated data rather than sessionwise?

% Let the data control which animal is being examined and then look over sessions for fit params within that.
all_sess_files = dir('*NAc*dat.mat');
figure(900); clf; unos = 1:3:34; dos = 2:3:35; tres = 3:3:36;

for zz=1:numel(all_sess_files)

    % load the corresponding hexa_data_an     
    T = load(all_sess_files(zz).name);

    % Examine the evolution of port choice over all sessions that a
    % given animals has
    gmap                = TNC_CreateRBColormap(8,'grima');

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
    ports               = T.hexa_data.port_n(uv_inds);

    rw_inds             = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==1);
    rtimes              = T.hexa_data.event_time_con(rw_inds);
    rports              = T.hexa_data.port_n(rw_inds);
    
% Get reward times with valid photometry
% also get a logical indexing from all visits/rewards to the subset for which I have photometry
    rw_p_inds           = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==1 & ~isnan(T.hexa_data.photo_i));
    rw_p_logic          = ismember(rw_inds,rw_p_inds);
    
    ur_p_inds           = find(T.hexa_data.unique_vis==1 & T.hexa_data.rewarded==0 & ~isnan(T.hexa_data.photo_i));
    rw_p_logic          = ismember(uv_inds,ur_p_inds);

    rw_p_inds_i         = T.hexa_data.photo_i(rw_p_inds);
    rw_p_inds_s         = T.hexa_data.session_n(rw_p_inds);
    rw_p_inds_p         = T.hexa_data.port_n(rw_p_inds);

    ur_p_inds_i         = T.hexa_data.photo_i(ur_p_inds);
    ur_p_inds_s         = T.hexa_data.session_n(ur_p_inds);
    ur_p_inds_p         = T.hexa_data.port_n(ur_p_inds);

% grab photometry events
    photo_event_win     = [240 400];
    photo_data          = sgolayfilt(T.hexa_data.photo.dFF,3,51);

    figure(701); clf; 
    sess_map = TNC_CreateRBColormap(6,'grima');
    
    for sess=unique(rw_p_inds_s)'
        figure(701);
        da_resp_rew(sess).da_sink= TNC_ExtTrigWins(photo_data(T.hexa_data.photo.sess==sess),rw_p_inds_i(rw_p_inds_s==sess),photo_event_win);
        da_resp_rew(sess).trap   = mean(da_resp_rew(sess).da_sink.wins(:,photo_event_win(1):photo_event_win(1)+120),2)-mean(da_resp_rew(sess).da_sink.wins(:,1:30),2);
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

            % subplot(1,5,5); plot(movmean(da_resp_rew(sess).trap,31),movmean(da_resp_ure(sess).trap,31),'color',[sess_map(sess+3,:) 1]/2,'LineWidth',2); hold on;

            figure(900);
            subplot(6,6,unos(zz));
            plot(da_resp_rew(sess).vist,movmean(da_resp_rew(sess).trap,31),'color',[sess_map(sess,:) 1],'LineWidth',2); hold on;
            plot(da_resp_ure(sess).vist,movmean(da_resp_ure(sess).trap,31),'color',[sess_map(sess,:) 1]/2,'LineWidth',2); hold on;
            xlabel('Session time (sec)'); xlim([0 180*5*60]);  box off; brk = strfind(all_sess_files(zz).name,'_'); title(all_sess_files(zz).name(1:brk-1));
            if zz==1
                ylabel('DA Response (Rew; Light) (UnRew; Dark)');
            end
    end


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

    % use this per_sess_rank and session_ids and visit_matrix to make concatenated reward schedule for entire animal's data
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
        plot(visit_matrix_sm(qqq,:),'color',gmap(qqq,:),'linewidth',2); hold on;
        xlim([0 size(visit_matrix,2)]); ylim([0 0.15]);
        ylabel('Fraction of visits');

        reward_matrix(qqq,ceil(rtimes(rports==qqq)))=1;
        reward_matrix_sm(qqq,:) = movmean(reward_matrix(qqq,:),30*60);

        figure(9); subplot(212);
        semilogy(reward_matrix_sm(qqq,:),'color',gmap(qqq,:),'linewidth',2); hold on;
        xlim([0 size(reward_matrix,2)]); ylim([0 0.075]);
        ylabel('Fraction of visits rewarded');

        figure(900);
        subplot(6,6,dos(zz));
        plot(visit_matrix_sm(qqq,:),'color',gmap(qqq,:),'linewidth',2); hold on;
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

    % alternative is to load optimal values for all available fit sessions
    % for this animal and concatenate
    dir_breaks = strfind(T.hexa_data.filename,'/');
    just_file = T.hexa_data.filename(dir_breaks(end)+1:end);
    mname_end = strfind(just_file,'_');
    mouse_name = just_file(1:mname_end-1);
    all_opt_files = dir([mouse_name '*opt*']);
    alpha = [];

    all_com = zeros(numel(all_opt_files),3);

    for zzz=1:numel(all_opt_files)

        S = load(all_opt_files(zzz).name);
            
        loss = S.opt_r2_tensor;
        
        [top_xperc_inds] = find(loss>frac*max(loss,[],"all"));

        [a2_inds,a4_inds,de_inds] = ind2sub(size(S.opt_r2_tensor),top_xperc_inds);
        [all_com(zzz,:)] = centerOfMass3D(S.a2_vec(a2_inds), S.a4_vec(a4_inds), S.de_vec(de_inds), loss(top_xperc_inds)');

        v_ind = 1:numel(find(sum(S.visit_matrix,1)==1));
        this_alpha = 0.001 + (all_com(zzz,1) ./ (1+exp((all_com(zzz,2)-v_ind)/(all_com(zzz,2)./6)))) .*  (all_com(zzz,1)*exp(-v_ind/(5*all_com(zzz,2))));

        alpha = [alpha this_alpha];
    end

    figure(12); clf; plot(times,alpha); xlim([0 max(times)]); box off;

    % doesn't really work at the moment - need to analyze each visit matrix independently
    num_iter            = 1; 
    vismat              = zeros(6,size(visit_matrix,2),num_iter);
    rewmat              = zeros(6,size(visit_matrix,2),num_iter);
    trans_r2_iter       = zeros(1,num_iter);
    income_r2_iter      = zeros(1,num_iter);
    income_model        = zeros(1,numel(income),num_iter);
    p_reward            = zeros(6,size(visit_matrix,2),num_iter);

    cost_per_port =                 ...
    [1	14	18	70	72.2	65.5;   ...
    14	1	22.8	56	65.5	42; ...
    18	22.8	1	72.2	70	56; ...
    70	56	72.2	1	18	22.8;   ...
    72.2	65.5	70	18	1	14; ...
    65.5	42	56	22.8	14	1];

    for iter = 1:num_iter
        [trans_r2_iter(iter), income_r2_iter(iter), vismat(:,:,iter), rewmat(:,:,iter), p_reward(:,:,iter), income_model(:,:,iter)] = HX_model_session_forAlphaConcat(alpha,visit_matrix,cost_per_port,rew_sched,income);
    end

    visits_for_LL = squeeze(mean(vismat,3));
        visits_for_LL_sm = zeros(size(visits_for_LL));
    rewards_for_LL = squeeze(mean(rewmat,3));
        rewards_for_LL_sm = zeros(size(rewards_for_LL));
    p_reward_all = squeeze(mean(p_reward,3));

    figure(8); clf;
    plot(1:numel(all_visits),[0 diff(all_sessid)],'k-'); hold on;
    plot(1:numel(all_visits),income,'color',gmap(3,:),'linewidth',2);
    plot(1:numel(all_visits),mean(income_model,3),'color',gmap(5,:),'linewidth',2);
    xlim([0 numel(all_visits)]); box off;

    figure(900);
    subplot(6,6,tres(zz));
    plot(1:numel(all_visits),[0 diff(all_sessid)],'k-'); hold on;
    plot(1:numel(all_visits),income,'color',gmap(3,:),'linewidth',2);
    plot(1:numel(all_visits),mean(income_model,3),'color',gmap(5,:),'linewidth',2);
    xlim([0 5.4e3]); box off;

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
        plot(visits_for_LL_sm(qqq,:),'color',gmap(qqq,:),'linewidth',2); hold on;
        xlim([0 size(visits_for_LL,2)]); ylim([0 0.15]);
        ylabel('Fraction of visits');

        rewards_for_LL(qqq,ceil(rtimes_m(rports_m==qqq)))=1;
        rewards_for_LL_sm(qqq,:) = movmean(rewards_for_LL(qqq,:),30*60);

        figure(11); subplot(212);
        semilogy(rewards_for_LL_sm(qqq,:),'color',gmap(qqq,:),'linewidth',2); hold on;
        xlim([0 size(rewards_for_LL,2)]); ylim([0 0.075]);
        ylabel('Fraction of visits rewarded');
    end
    figure(11); subplot(211); yyaxis right; plot(session_ids,'k'); box off;
    ylabel('Session');
    all_rewards = sum(rewards_for_LL_sm,1);
    figure(11); subplot(212); plot(all_rewards,'k-'); yyaxis right; plot(session_ids,'k'); box off;
    ylabel('Session');
    
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

    % figure(100); clf;
    % histogram(null_predictions);
    % hold on;
    % text(choice_predictions,5,'v');

    sh_nom_pvalue = numel(find(null_predictions<choice_predictions)) ./ numel(null_predictions)

    % align to session 2->3 transition
    transition = find(session_ids==2 & [diff(session_ids) 0]==1);
    if numel(transition)>0
        twin       = [ 2000 5000 ];
        figure(101); clf;
        subplot(121);
        plot([0 0 0 0.08],'k-','linewidth',2); hold on;
        ylabel('P(choice)'); xlabel('Time from switch (s)'); title('Model');
        subplot(122);
        plot([0 0 0 0.08],'k-','linewidth',2); hold on;
        ylabel('P(choice)'); xlabel('Time from switch (s)'); title('Data');
        for qqq = 1:6
            subplot(121);
            plot(-twin(1):twin(2),visits_for_LL_sm(qqq,transition-twin(1):transition+twin(2)),'color',gmap(qqq,:),'linewidth',2); hold on;
            box off;
            subplot(122);
            plot(-twin(1):twin(2),visit_matrix_sm(qqq,transition-twin(1):transition+twin(2)),'color',gmap(qqq,:),'linewidth',2); hold on;
            box off;
        end
    end

end