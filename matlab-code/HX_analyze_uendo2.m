% HX_analyze_microendo

%% Data Loading / Construction of data matrix
clear trace_mat

filename_traces         = '6PM4_2022-10-14_celltraces.csv';
filename_traces_props   = '6PM4_2022-10-14_celltraces-props.csv';
filename_behav          = '6PM4_behaviour.csv';

exag    = TNC_CreateRBColormap(1000,'exag');
catmap  = TNC_CreateRBColormap(1000,'cat2');

traces          = readtable(filename_traces);
traces_props    = readtable(filename_traces_props);
behav           = readtable(filename_behav);

trace_mat.dff   = [];

for kk=1:numel(traces.Properties.VariableNames)
    
    if numel(strfind(traces.Properties.VariableNames{kk},'C'))>0
        tmp = traces.(traces.Properties.VariableNames{kk});
        trace_mat.dff = [trace_mat.dff ; tmp'];
    elseif numel(strfind(traces.Properties.VariableNames{kk},'Var1'))
        trace_mat.t = traces.(traces.Properties.VariableNames{kk});
    end
    
end

trace_mat.x = traces_props.CentroidX;
trace_mat.y = traces_props.CentroidY;

mappedX = tsne(trace_mat.dff);

[idx,C,sumdist] = kmeans(mappedX,5,'Distance','cityblock', 'Display','final','Replicates',10);
[~,inds]=sort(idx);

figure(1); clf;
    imagesc(trace_mat.dff(inds,:),[0 10]); colormap(exag);

figure(2); clf;
    scatter(mappedX(:,1),mappedX(:,2),100,idx,'filled'); colormap(catmap);
    axis square;

figure(3); clf;
    imagesc(corr(trace_mat.dff(inds,:)'),[0 1]); colormap(exag);
    axis square;

figure(4); clf;
    scatter(trace_mat.x,trace_mat.y,50,idx,'filled'); colormap(catmap);
    axis([0 600 0 600]);
    axis square;


%% Load behavioral tracking data from raw files

filename_track = '../raw/6PM4_hex_foraging_V2_2022-10-14-102642.csv';
tracker = readtable(filename_track);

track_dat.pos.y = tracker.YLocation .* 0.3;
fill = find( track_dat.pos.y == 0 );
cont = find( track_dat.pos.y ~= 0 );

track_dat.pos.y(fill) = interp1(cont,track_dat.pos.y(cont),fill);

tmp_img_x = tracker.XLocation.* 0.3;
tmp_img_x(fill) = interp1(cont,tmp_img_x(cont),fill);

track_dat.pos.x = tracker.StagePostion.* 0.0069 + tmp_img_x;

figure(101); clf;
plot(sgolayfilt( track_dat.pos.x , 3 , 101) , sgolayfilt( track_dat.pos.y , 3, 101) , 'color', [0 0 0 0.1]);

debug = 0;
filename_gpio = '../raw/6PM4_2022-10-14_GPIO.csv';
gpio = readtable(filename_gpio);

tf = contains(gpio.ChannelName,'GPIO-1');
track_dat.val = gpio.Value(tf==1)';
times = gpio.Time_s_(tf==1)';

pulses = find(track_dat.val>1e4 & [0 diff(track_dat.val)]>0);
start_pulses = pulses(find([251 diff(pulses)>250]));


if debug

    figure(100); clf;
    subplot(141);
    plot(tracker.StagePostion.* 0.0069); 
    yyaxis right;
    plot(tracker.XLocation.* 0.3);
    subplot(142);
    plot(tracker.YLocation.* 0.3);
    subplot(143);
    plot(sgolayfilt( tracker.StagePostion.* 0.0069 + tracker.XLocation.* 0.3 , 3 , 101));
    subplot(144);
    plot(sgolayfilt( behav.pos.x , 3 , 101) , sgolayfilt( behav.pos.y , 3, 101) );
    
    figure(8); clf; plot(track_dat.val);
    hold on;
    plot(pulses,ones(1,numel(pulses)).*1e4,'r*');
    plot(start_pulses,ones(1,numel(start_pulses)).*1e4,'go');

end

track_dat.gpio.times = times(start_pulses);

%% Get temporal alignment
% can be aligned by comparing to tracker.Ext_ClockCount which increments by
% 1 each time a GPIO pulse is received

figure(1); clf; 
subplot(1,6,1);
plot(track_dat.gpio.times, find(diff(tracker.Ext_ClockCount)==1), 'k-');

% fit this curve and use offset and slope to compute alignment
t_off_fit = polyfit(track_dat.gpio.times, find(diff(tracker.Ext_ClockCount)==1), 1);

hold on; 
% plot(tracker.FrameNumber,'k-'); hold on;
plot(trace_mat.t,polyval(t_off_fit,trace_mat.t),'ro');

% this fit should provide me with a way to generate the appropriate sample
% from position for any given index in the imaging data - right? double
% check this...
track_dat.pos.xs = sgolayfilt( track_dat.pos.x , 3 , 51);
track_dat.pos.ys = sgolayfilt( track_dat.pos.y , 3 , 51);

track_dat.vis.uv = zeros(1,numel(track_dat.pos.x));
track_dat.vis.rw = zeros(1,numel(track_dat.pos.x));

track_dat.vis.uv(unique(behav.vid_i(behav.unique_visit==1))) = 1;
track_dat.vis.rw(unique(behav.vid_i(behav.rewarded==1))) = 1;

% id ports from visit times
figure(50); clf;
subplot(1,7,1:5);
plot(track_dat.pos.xs,track_dat.pos.ys,'Color',[0 0 0 0.1]); hold on;
scatter(track_dat.pos.x(track_dat.vis.uv==1),track_dat.pos.y(track_dat.vis.uv==1),25,'k','filled')
scatter(track_dat.pos.x(track_dat.vis.rw==1),track_dat.pos.y(track_dat.vis.rw==1),25,'r')
axis off;

port_cmap = TNC_CreateRBColormap(8,'grima');
for pp = 1:6
    subplot(1,7,6);
    plot(cumsum(trace_mat.pid==pp),'color',port_cmap(pp,:),'LineWidth',2); hold on;
    subplot(1,7,7);
    tmp = trace_mat.rew;
    tmp(trace_mat.pid~=pp) = 0;
    plot(cumsum(tmp),'color',port_cmap(pp,:),'LineWidth',2); hold on;    
end
subplot(1,7,6); ylabel('\Sigma Visits'); legend({'30' '60' '240' '240' '1200' '2400'},'Location','northwest');
subplot(1,7,7); ylabel('\Sigma Rewards');

frames = polyval(t_off_fit,trace_mat.t);
valid_frames = find(frames>0.5 & frames<=tracker.FrameNumber(end));

positions_per_imageFrame.x = track_dat.pos.xs( round(frames(valid_frames)) );
positions_per_imageFrame.y = track_dat.pos.ys( round(frames(valid_frames)) );

% for visits need to find the closest frame to each unique_visit
trace_mat.uv    = zeros(1,numel(trace_mat.t));
trace_mat.rew   = zeros(1,numel(trace_mat.t));
trace_mat.pid   = zeros(1,numel(trace_mat.t));
all_uv = find(behav.unique_visit==1);
for jj=all_uv'
    this_time = behav.micro_i(jj);
    this_rew = behav.rewarded(jj);
    this_pid = behav.port(jj);
    frame_i = find(trace_mat.t>=this_time,1,'first');
    trace_mat.uv(frame_i)=1;
    trace_mat.rew(frame_i)=this_rew;
    trace_mat.pid(frame_i)=this_pid;
end

subplot(1,6,2:6);
plot(track_dat.pos.xs,track_dat.pos.ys,'k-');
hold on;
plot(positions_per_imageFrame.x,positions_per_imageFrame.y,'r.');

%% try a 3d plot of activity vs position for quick visualization

figure(10); clf;

subplot(121);
plot3(positions_per_imageFrame.x,positions_per_imageFrame.y,valid_frames,'k','Color',[0 0 0 0.2]);
hold on;
plot3(positions_per_imageFrame.x(trace_mat.dff(inds(91),valid_frames)>2),positions_per_imageFrame.y(trace_mat.dff(inds(91),valid_frames)>2),valid_frames(trace_mat.dff(inds(91),valid_frames)>2),'r.');

zlabel('Time'); xlabel('X pos'); ylabel('Y pos');

subplot(122);
plot3(positions_per_imageFrame.x,positions_per_imageFrame.y,valid_frames,'k','Color',[0 0 0 0.2]);
hold on;
plot3(positions_per_imageFrame.x(trace_mat.dff(inds(140),valid_frames)>2),positions_per_imageFrame.y(trace_mat.dff(inds(140),valid_frames)>2),valid_frames(trace_mat.dff(inds(140),valid_frames)>2),'b.');

zlabel('Time'); xlabel('X pos'); ylabel('Y pos');

%% Segment into trajectories relative to position of rank1 port

[~,positions_per_imageFrame.d] = cart2pol( ...
    positions_per_imageFrame.x-mode(positions_per_imageFrame.x), ...
    positions_per_imageFrame.y-430);


figure(11); clf;
plot(positions_per_imageFrame.d);

events = TNC_ExtractPeaks(sgolayfilt(positions_per_imageFrame.d',3,21),150,75,1);

% [sink] = TNC_ExtTrigWins3d(trace_mat.dff,events.stops,[500 150]);
% figure(2); clf; imagesc(sink.win_avg(perm,:));

[sink.x] = TNC_ExtTrigWins(positions_per_imageFrame.x,events.stops,[500 150]);
[sink.y] = TNC_ExtTrigWins(positions_per_imageFrame.y,events.stops,[500 150]);

delta_y = mean(sink.y.wins(:,100:450),2)-450;
[~,sorty] = sort(delta_y);

figure(3); clf; subplot(131); imagesc(sink.x.wins(sorty,:)); 
subplot(132); imagesc(sink.y.wins(sorty,:)); 

[sink] = TNC_ExtTrigWins3d(trace_mat.dff,find(trace_mat.rew==1),[50 150]);
% figure(2); clf; imagesc(sink.win_avg(perm,:));

pid_rew = trace_mat.pid(find(trace_mat.rew==1));
[ps,sortp] = sort(pid_rew);

figure(4); clf;
for zz=1:144
    figure(4); subplot(12,12,zz);
    imagesc(squeeze(sink.wins(perm(zz),:,sortp))'); 
    colormap(flipud(bone)); axis off;
end

figure(20); clf;
subplot(1,8,1:5)
plot(trace_mat.uv.*75,'k','color',[0.7 0.7 0.7]);
hold on;
plot(trace_mat.rew.*75,'k');
plot(trace_mat.dff(perm(26),:));
plot(trace_mat.dff(perm(110),:));
plot(trace_mat.dff(perm(94),:));

subplot(186)
imagesc(squeeze(sink.wins(perm(26),:,sortp))'); 
colormap(flipud(bone)); ylabel('sorted by port id');
subplot(187)
imagesc(squeeze(sink.wins(perm(110),:,sortp))'); 
colormap(flipud(bone)); ylabel('sorted by port id');
subplot(188)
imagesc(squeeze(sink.wins(perm(94),:,sortp))'); 
colormap(flipud(bone)); ylabel('sorted by port id');


%% Try to look for place fields in a semistandard way

% Unwrap the 2D environment into a 1D array and align over time
xbins = 0:25:2200;
ybins = -25:25:525;
[positions_per_imageFrame.xd] = discretize(positions_per_imageFrame.x,xbins);
[positions_per_imageFrame.yd] = discretize(positions_per_imageFrame.y,ybins);

place_cell_map = TNC_CreateRBColormap(1024,'grima-denser');

trace_mat.dffxPos = trace_mat.dff(:,valid_frames);

mean_resp = zeros(numel(ybins),numel(xbins),size(trace_mat.dffxPos,1));

for pp = unique(positions_per_imageFrame.xd)'
    for qq = unique(positions_per_imageFrame.yd)'
        for zz = 1:numel(inds)
            this_loc_inds = find(positions_per_imageFrame.xd==pp & positions_per_imageFrame.yd==qq);
            mean_resp(qq,pp,zz) = mean(trace_mat.dffxPos(inds(zz),this_loc_inds),'omitnan');
        end
    end
end


pc_kern = TNC_CreateGaussian(100,1.5,200,1);
smth_resp = mean_resp;
figure(102); clf;
for zz = 1:numel(inds)
    subplot(13,12,zz);
    tmp = mean_resp(:,:,perm(zz));
    tmp(isnan(tmp)) = 0;
    % imagesc(xbins,ybins,conv2(tmp,pc_kern'*pc_kern,'same'),[0 max(conv2(tmp,pc_kern'*pc_kern,'same'),[],'all')]); colormap(place_cell_map); axis off; 
    % title(num2str(max(conv2(tmp,pc_kern'*pc_kern,'same'),[],'all')))
    smth_resp(:,:,perm(zz)) = conv2(tmp,pc_kern'*pc_kern,'same');
    mags(zz) = max(smth_resp(:,:,perm(zz)),[],'all');
    drawnow;
end

%% just the strongest responses
pc_kern = TNC_CreateGaussian(100,1.5,200,1);
figure(103); clf;
cnt=1;
for zz = find(mags>2.5)
    subplot(8,5,cnt);
    imagesc(xbins,ybins,smth_resp(:,:,perm(zz)),[0 mags(zz)]); 
    colormap(place_cell_map); axis off; 
    drawnow;
    cnt = cnt+1;
end

    %% Try to find sequences of activity

    ord_map = TNC_CreateRBColormap(146,'cpb');

    % first compute the correlation matrix 
    C = corr(trace_mat.dff');

    % Convert to a distance
    D = 1-abs(C);    
    D_vec = squareform(D,'tovector');
    Z = linkage(D_vec,'average');
    leafOrder = optimalleaforder(Z,D_vec);
    perm = leafOrder;

    C_ordered = C(perm, perm);

    figure(1); clf; 
    subplot(2,6,1); imagesc(C);
    subplot(2,6,2:6); imagesc(trace_mat.dff,[0 10]); xlim([0 1e4]);

    subplot(267); imagesc(C_ordered);
    subplot(2,6,8:12); imagesc(trace_mat.dff(perm,:),[0 10]); colormap(flipud(bone)); xlim([0 1e4]);

    % an interesting cluster is 
    clust_inds = perm %perm(85:135)
    range = 2.7e4:5.75e4;

    figure(2); clf;
    for vv=clust_inds

        plot(range,trace_mat.dff(vv,range)+(10*find(vv==clust_inds)),'color',ord_map(find(vv==clust_inds),:)); hold on; xlim([min(range) max(range)]);

    end
    plot(range,trace_mat.uv(range).*1500,'k-','Color',[0.75 0.75 0.75]);
    plot(range,trace_mat.rew(range).*1600,'k-');
    plot(range,50.*trace_mat.pid(range)+1600,'b-');

    figure(3); clf;
    plot(valid_frames,positions_per_imageFrame.x,'k','Color',[0 0 0 1]); xlim([0 5e4]);
    hold on;
    plot(valid_frames,positions_per_imageFrame.y,'k','Color',[1 0 0 1]); xlim([0 5e4]);
    plot(valid_frames,5*sum(trace_mat.dffxPos(perm,:),1),'k','Color',[0 0.67 1 1]); axis tight; xlim([0 5e4]); 
    