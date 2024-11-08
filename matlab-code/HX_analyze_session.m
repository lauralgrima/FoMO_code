function [hexa_data_an] = HX_analyze_session(hexa_data,session,photo_flag)

all_rew_inds            = find(hexa_data.rewarded==1 & ismember(hexa_data.session_n,session));
all_rew_ports           = hexa_data.port_n(all_rew_inds);

all_vis_inds            = find(hexa_data.unique_vis==1 & ismember(hexa_data.session_n,session));
all_vis_ports           = hexa_data.port_n(all_vis_inds);

all_vis_event_t         = hexa_data.event_time(all_vis_inds);
all_rew_event_t         = hexa_data.event_time(all_rew_inds);

hexa_data_an.photo_i    = hexa_data.photo_i(all_vis_inds);
hexa_data_an.video_i    = hexa_data.video_i(all_vis_inds);

hexa_data_an.filename   = hexa_data.filename;
hexa_data_an.session    = session;

figure(57); clf;
cat_map = TNC_CreateRBColormap(8,'mapb');

if numel(session)>1
    hexa_data_an.visits     = zeros(6,ceil(max(hexa_data.event_time_con(find(ismember(hexa_data.session_n,session))))));
    hexa_data_an.rewards    = zeros(6,ceil(max(hexa_data.event_time_con(find(ismember(hexa_data.session_n,session))))));
    hexa_data_an.sessID     = ones(1,ceil(max(hexa_data.event_time_con(find(ismember(hexa_data.session_n,session))))));
    for zz=session
        this_sess_start = ceil(max(hexa_data.event_time_con(find(hexa_data.session_n==zz,1,'first'))));
        this_sess_end = ceil(max(hexa_data.event_time_con(find(ismember(hexa_data.session_n,zz)))));
        hexa_data_an.sessID(this_sess_start:this_sess_end) = zz;
    end
else
    hexa_data_an.visits     = zeros(6,ceil(max(hexa_data.event_time(find(ismember(hexa_data.session_n,session))))));
    hexa_data_an.rewards    = zeros(6,ceil(max(hexa_data.event_time(find(ismember(hexa_data.session_n,session))))));
    hexa_data_an.sessID     = session * ones(1,ceil(max(hexa_data.event_time_con(find(ismember(hexa_data.session_n,session))))));
end
hexa_data_an.vi.avg = zeros(1,6);
hexa_data_an.vi.std = zeros(1,6);
hexa_data_an.rw.avg = zeros(1,6);
hexa_data_an.rw.std = zeros(1,6);

for qq=unique(all_vis_ports)'

    hexa_data_an.port_rew(qq).ts = hexa_data.event_time(all_rew_inds(all_rew_ports==qq));
    hexa_data_an.port_rew(qq).iti = diff(hexa_data_an.port_rew(qq).ts);
    disp(['Port ' num2str(qq) ' mean ri: ' num2str(mean(hexa_data_an.port_rew(qq).iti)) ' +/- ' num2str(std(hexa_data_an.port_rew(qq).iti))]);

    subplot(131);
    hexa_data_an.port_rew(qq).histX = 0.2:0.2:4.4;
    hexa_data_an.port_rew(qq).histY = hist(log10(hexa_data_an.port_rew(qq).iti),hexa_data_an.port_rew(qq).histX);
    plot(hexa_data_an.port_rew(qq).histX , cumsum(hexa_data_an.port_rew(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','2','3','4','5','6'},'Location','northwest');
        xlabel('Log10 Interreward Interval'); ylabel('Count'); box off;
    end

    hexa_data_an.rewards(qq,ceil(hexa_data.event_time(all_rew_inds(all_rew_ports==qq)))) = 1;
    hexa_data_an.rw.avg(qq) = mean(hexa_data_an.port_rew(qq).iti);
    hexa_data_an.rw.std(qq) = std(hexa_data_an.port_rew(qq).iti);

    hexa_data_an.port_vis(qq).ts = hexa_data.event_time(all_vis_inds(all_vis_ports==qq));
    hexa_data_an.port_vis(qq).iti = diff(hexa_data_an.port_vis(qq).ts);
    disp(['Port ' num2str(qq) ' mean vi: ' num2str(mean(hexa_data_an.port_vis(qq).iti)) ' +/- ' num2str(std(hexa_data_an.port_vis(qq).iti))]);

    subplot(132)
    hexa_data_an.port_vis(qq).histX = 0.2:0.2:4.4;
    hexa_data_an.port_vis(qq).histY = hist(log10(hexa_data_an.port_vis(qq).iti),hexa_data_an.port_vis(qq).histX);
    plot(hexa_data_an.port_vis(qq).histX , cumsum(hexa_data_an.port_vis(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','2','3','4','5','6'},'Location','northwest');
        xlabel('Log10 Intervisit Interval'); ylabel('Count'); box off;
    end
    
    hexa_data_an.visits(qq,ceil(hexa_data.event_time(all_vis_inds(all_vis_ports==qq)))) = 1;
    hexa_data_an.vi.avg(qq) = mean(hexa_data_an.port_vis(qq).iti);
    hexa_data_an.vi.std(qq) = std(hexa_data_an.port_vis(qq).iti);

end

hexa_data_an.all_rew_inds=all_rew_inds;
hexa_data_an.all_rew_ports=all_rew_ports;
hexa_data_an.all_vis_inds=all_vis_inds;
hexa_data_an.all_vis_ports=all_vis_ports;

subplot(133);
plot([1 max(log10(hexa_data_an.rw.avg))+1],[1 max(log10(hexa_data_an.rw.avg))+1],'k--'); hold on;
scatter(log10(hexa_data_an.vi.avg),log10(hexa_data_an.rw.avg),50,'filled');
axis([1 max(log10(hexa_data_an.rw.avg))+1 1 max(log10(hexa_data_an.rw.avg))+1]);
ylabel('Inter Reward Interval');
xlabel('Inter Visit Interval');
box off;

%%
smth_type = 1; % 1: symmetric gaussian; 2: history gaussian; 3: full history
highpass_photo = 1;
trial_win = 25;
trial_kernel = TNC_CreateGaussian(500,trial_win,1000,1);
if smth_type == 2
    trial_kernel(500:end)=0;
end

if photo_flag
    % get photometry events for every visit with corresponding vectors
    % reflecting which port and whether rewarded
    port_color_map      = TNC_CreateRBColormap(8,'mapb');
    sym_map             = TNC_CreateRBColormap(1000,'bb-sym');    
    photo_event_win     = [250 500];

    if numel(session)>1
        photo_sess_data     = sgolayfilt(hexa_data.photo.dFF(ismember(hexa_data.photo.sess,session)),3,51);
        visit_indices       = find(hexa_data.unique_vis==1 & ~isnan(hexa_data.photo_i_con) & ismember(hexa_data.session_n,session));
        photo_visit_inds    = hexa_data.photo_i_con(visit_indices);
    else
        photo_sess_data     = sgolayfilt(hexa_data.photo.dFF(hexa_data.photo.sess==session),3,51);
        visit_indices       = find(hexa_data.unique_vis==1 & ~isnan(hexa_data.photo_i) & hexa_data.session_n==session);
        photo_visit_inds    = hexa_data.photo_i(visit_indices);
    end
    % need a correction for the fact that photo_i resets each session
    % whereas other indexing operations do not
    
    port_visit_ids      = hexa_data.port_n(visit_indices);
    rew_visit_ids       = hexa_data.rewarded(visit_indices);
    rew_rate_in_visits  = conv(rew_visit_ids,[0 ones(1,50) 0]/50,'same');

    photo_visit_inds(1)
    photo_visit_inds(end)
    numel(photo_sess_data)

    if highpass_photo
        d2 = designfilt("highpassiir",FilterOrder=12, ...
            HalfPowerFrequency=0.0001,DesignMethod="butter");
        yH = filtfilt(d2,photo_sess_data);
        photo_sess_data = yH;
    end


    sink = TNC_ExtTrigWins(photo_sess_data,photo_visit_inds,photo_event_win);
    [~,sort_inds] = sort(port_visit_ids.*(rew_visit_ids-0.5),'descend');

    figure(301); clf;
    imagesc(sink.wins(sort_inds,:),[-5 5]); colormap(sym_map);

    hexa_data_an.da_visit_resp  = sink;
    hexa_data_an.da_visit_ids   = port_visit_ids;
    hexa_data_an.da_visit_rew   = rew_visit_ids;

    hexa_data_an.da_hand1 = figure(300); clf;
    for pp=1:6
        subplot(1,6,pp);

        pN_inds = find( port_visit_ids==pp & rew_visit_ids==0 );
        shadedErrorBar(sink.range,mean(sink.wins(pN_inds,:),1),std(sink.wins(pN_inds,:),[],1)./sqrt(numel(pN_inds)),{'color',[0.5 0.5 0.5]}); hold on;
        
        pR_inds = find( port_visit_ids==pp & rew_visit_ids==1 );        
        shadedErrorBar(sink.range,mean(sink.wins(pR_inds,:),1),std(sink.wins(pR_inds,:),[],1)./sqrt(numel(pR_inds)),{'color',port_color_map(pp,:)}); hold on;
        axis([sink.range(1) sink.range(end) -2 4]); box off; grid on;
    end

    hexa_data_an.da_hand2 = figure(302); clf;
    if numel(session)>1
        for pp=1:6
            p_rew_this_port = cumsum(rew_visit_ids.*port_visit_ids==pp) ./ cumsum(port_visit_ids==pp);
    
            subplot(1,6,pp);
            plot([0 max(hexa_data.event_time_con(visit_indices))],[0 0],'k-');  hold on;
    
            pN_inds = find( port_visit_ids==pp & rew_visit_ids==0 );
            scatter(hexa_data.event_time_con(visit_indices(pN_inds)),min(sink.wins(pN_inds,photo_event_win(1):photo_event_win(1)+300),[],2),50,port_color_map(pp,:),'MarkerEdgeAlpha',0.2); colormap(port_color_map);
            plot(hexa_data.event_time_con(visit_indices(pN_inds)),conv( min(sink.wins(pN_inds,photo_event_win(1):photo_event_win(1)+300),[],2) , [0 ones(1,5) 0]/5,'same'),'color',port_color_map(pp,:),'LineWidth',1); colormap(port_color_map);
    
            pR_inds = find( port_visit_ids==pp & rew_visit_ids==1 );
    
            scatter(hexa_data.event_time_con(visit_indices(pR_inds)),max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+300),[],2),50,port_color_map(pp,:),"filled",'MarkerFaceAlpha',0.2); colormap(port_color_map);
            if numel(pR_inds)>5
                plot(hexa_data.event_time_con(visit_indices(pR_inds)),conv( max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+300),[],2) , [0 ones(1,5) 0]/5,'same'),'color',port_color_map(pp,:),'LineWidth',2); colormap(port_color_map);
            else
                plot(hexa_data.event_time_con(visit_indices(pR_inds)),max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+325),[],2),'color',port_color_map(pp,:),'LineWidth',2); colormap(port_color_map);
            end
            axis([0 max(hexa_data.event_time_con(visit_indices)) -3 6]);        
    
            box off; grid on; title(['Port n=' num2str(pp)]);
            
            if pp==1
                ylabel('Integrated DA resp.'); xlabel('Unique visits'); 
            end
            yyaxis right;
            plot(hexa_data.event_time_con(visit_indices),p_rew_this_port','color',[0 0 0 0.2],'LineWidth',2); 
            axis([0 max(hexa_data.event_time_con(visit_indices)) 0 1]);
            if pp==1
                ylabel('P(rew|visit)');
            end
        end
    else
        for pp=1:6
            p_rew_this_port = cumsum(rew_visit_ids.*port_visit_ids==pp) ./ cumsum(port_visit_ids==pp);
    
            subplot(1,6,pp);
            plot([0 max(hexa_data.event_time(visit_indices))],[0 0],'k-');  hold on;
    
            pN_inds = find( port_visit_ids==pp & rew_visit_ids==0 );
            scatter(hexa_data.event_time(visit_indices(pN_inds)),min(sink.wins(pN_inds,photo_event_win(1):photo_event_win(1)+300),[],2),50,port_color_map(pp,:),'MarkerEdgeAlpha',0.2); colormap(port_color_map);
            plot(hexa_data.event_time(visit_indices(pN_inds)),conv( min(sink.wins(pN_inds,photo_event_win(1):photo_event_win(1)+300),[],2) , [0 ones(1,5) 0]/5,'same'),'color',port_color_map(pp,:),'LineWidth',1); colormap(port_color_map);
    
            pR_inds = find( port_visit_ids==pp & rew_visit_ids==1 );
    
            scatter(hexa_data.event_time(visit_indices(pR_inds)),max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+300),[],2),50,port_color_map(pp,:),"filled",'MarkerFaceAlpha',0.2); colormap(port_color_map);
            if numel(pR_inds)>5
                plot(hexa_data.event_time(visit_indices(pR_inds)),conv( max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+300),[],2) , [0 ones(1,5) 0]/5,'same'),'color',port_color_map(pp,:),'LineWidth',2); colormap(port_color_map);
            else
                plot(hexa_data.event_time(visit_indices(pR_inds)),max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+300),[],2),'color',port_color_map(pp,:),'LineWidth',2); colormap(port_color_map);
            end
            axis([0 max(hexa_data.event_time(visit_indices)) -3 6]);        
    
            box off; grid on; title(['Port n=' num2str(pp)]);
            
            if pp==1
                ylabel('Integrated DA resp.'); xlabel('Unique visits'); 
            end
            yyaxis right;
            plot(hexa_data.event_time(visit_indices),p_rew_this_port','color',[0 0 0 0.2],'LineWidth',2); 
            axis([0 max(hexa_data.event_time(visit_indices)) 0 1]);
            if pp==1
                ylabel('P(rew|visit)');
            end
        end
    end

    figure(299); clf; clear *_all;

    hexa_data_an.visit_indices = visit_indices;
    da_kernel = TNC_CreateGaussian(500,2,1000,1);

    for pp=1:6

        if smth_type==3
            p_rew_all(pp,:) = cumsum(rew_visit_ids.*port_visit_ids==pp) ./ cumsum(port_visit_ids==pp);
            p_choice_all(pp,:) = cumsum(port_visit_ids==pp) ./ cumsum(port_visit_ids>0);
        else
            p_rew_all(pp,:) = conv( (rew_visit_ids.*port_visit_ids==pp) , trial_kernel , 'same' );
            p_choice_all(pp,:) = conv( (port_visit_ids==pp) , trial_kernel , 'same' );
        end

        pR_inds = find( port_visit_ids==pp & rew_visit_ids==1 );
        da_resp(pp).int = trapz(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+325),2);
        da_resp(pp).t   = hexa_data.event_time_con(visit_indices(pR_inds));
    end

    hexa_data_an.da_resp = da_resp;
    hexa_data_an.da_resp_data = sink;
    hexa_data_an.da_resp_data.p_vis_id = port_visit_ids;
    hexa_data_an.da_resp_data.pp_vis_id = [1 ; port_visit_ids(1:end-1)];
    hexa_data_an.da_resp_data.r_vis_id = rew_visit_ids;

    pN_inds = find( port_visit_ids>0 & rew_visit_ids==0 );
    hexa_data_an.da_resp_all.u  = min(sink.wins(pN_inds,photo_event_win(1):photo_event_win(1)+350),[],2)-mean(sink.wins(pN_inds,1:photo_event_win(1)-100),2);
    hexa_data_an.da_resp_all.tu = hexa_data.event_time(visit_indices(pN_inds));
    hexa_data_an.da_resp_all.iu = pN_inds;
    hexa_data_an.da_resp_all.pu = port_visit_ids(pN_inds);
    hexa_data_an.da_resp_all.ppu= [0 ; port_visit_ids(pN_inds(2:end)-1)];
    
    pR_inds = find( port_visit_ids>0 & rew_visit_ids==1 );
    hexa_data_an.da_resp_all.r  = max(sink.wins(pR_inds,photo_event_win(1):photo_event_win(1)+350),[],2)-mean(sink.wins(pR_inds,1:photo_event_win(1)-100),2);
    hexa_data_an.da_resp_all.t  = hexa_data.event_time(visit_indices(pR_inds));
    hexa_data_an.da_resp_all.i  = pR_inds;
    hexa_data_an.da_resp_all.p  = port_visit_ids(pR_inds);
    hexa_data_an.da_resp_all.pp = port_visit_ids(pR_inds(pR_inds>1)-1);
    hexa_data_an.da_resp_all.tl = [ 1e3 diff(hexa_data.event_time(visit_indices(pR_inds))') ];

    if numel(session)>1
        subplot(131);
        plot(hexa_data.event_time_con(visit_indices),mean(p_choice_all,1),'linewidth',3,'color',[0 0 0]); hold on;
        for pp=1:6
            plot(hexa_data.event_time_con(visit_indices),p_choice_all(pp,:),'linewidth',2,'color',port_color_map(pp,:)); 
        end
        axis([0 max(hexa_data.event_time_con(visit_indices)) 0 0.5]);
        ylabel(['P(visit,port) | sigma=' num2str(trial_win)]); box off;
        hexa_data_an.p_choice_all = p_choice_all;
    
        subplot(132);
        plot(hexa_data.event_time_con(visit_indices),mean(p_rew_all,1),'linewidth',3,'color',[0 0 0]); hold on;
        for pp=1:6
            plot(hexa_data.event_time_con(visit_indices),p_rew_all(pp,:),'linewidth',2,'color',port_color_map(pp,:)); 
        end
        axis([0 max(hexa_data.event_time_con(visit_indices)) 0 0.35]);
        ylabel(['P(rew,port) | sigma=' num2str(trial_win)]); box off;
        hexa_data_an.p_rew_all = p_rew_all;
        
        subplot(133);
        for pp=1:6
            plot(da_resp(pp).t,conv(da_resp(pp).int,da_kernel,'same'),'linewidth',2,'color',port_color_map(pp,:)); hold on;
        end
        axis([0 max(hexa_data.event_time_con(visit_indices)) -250 1000]);
        ylabel('DA resp.'); box off;
        yyaxis right;
        plot(hexa_data.event_time_con(visit_indices),conv( [1 diff(hexa_data.event_time_con(visit_indices))']  , trial_kernel , 'same' ) ,'color' , [0 0 0 0.1], 'linewidth' , 4);
        
        [~,rank1_port] = max(mean(p_rew_all,2)); 
        pR_inds = find( port_visit_ids==rank1_port & rew_visit_ids==1);

        hexa_data_an.event_time = hexa_data.event_time_con;
        
    else
        subplot(131);
        plot(hexa_data.event_time(visit_indices),mean(p_choice_all,1),'linewidth',3,'color',[0 0 0]); hold on;
        for pp=1:6
            plot(hexa_data.event_time(visit_indices),p_choice_all(pp,:),'linewidth',2,'color',port_color_map(pp,:)); 
        end
        axis([0 max(hexa_data.event_time(visit_indices)) 0 0.5]);
        ylabel(['P(visit,port) | sigma=' num2str(trial_win)]); box off;
        hexa_data_an.p_choice_all = p_choice_all;
    
        subplot(132);
        plot(hexa_data.event_time(visit_indices),mean(p_rew_all,1),'linewidth',3,'color',[0 0 0]); hold on;
        for pp=1:6
            plot(hexa_data.event_time(visit_indices),p_rew_all(pp,:),'linewidth',2,'color',port_color_map(pp,:)); 
        end
        axis([0 max(hexa_data.event_time(visit_indices)) 0 0.35]);
        ylabel(['P(rew,port) | sigma=' num2str(trial_win)]); box off;
        hexa_data_an.p_rew_all = p_rew_all;
        
        subplot(133);
        for pp=1:6
            plot(da_resp(pp).t,conv(da_resp(pp).int,da_kernel,'same'),'linewidth',2,'color',port_color_map(pp,:)); hold on;
        end
        axis([0 max(hexa_data.event_time(visit_indices)) -250 1000]);
        ylabel('DA resp.'); box off;
        yyaxis right;
        plot(hexa_data.event_time(visit_indices),conv( [1 diff(hexa_data.event_time(visit_indices))']  , trial_kernel , 'same' ) ,'color' , [0 0 0 0.1], 'linewidth' , 4);
        
        [~,rank1_port] = max(mean(p_rew_all,2)); 
        pR_inds = find( port_visit_ids==rank1_port & rew_visit_ids==1);

        hexa_data_an.event_time = hexa_data.event_time;

    end
end




