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

frames = polyval(t_off_fit,trace_mat.t);
valid_frames = find(frames>0.5 & frames<=tracker.FrameNumber(end));

positions_per_imageFrame.x = track_dat.pos.xs( round(frames(valid_frames)) );
positions_per_imageFrame.y = track_dat.pos.ys( round(frames(valid_frames)) );

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

%% Try to look for place fields in a semistandard way

% Unwrap the 2D environment into a 1D array and align over time
xbins = 50:10:2200;
ybins = 0:10:500;
[positions_per_imageFrame.xd] = discretize(positions_per_imageFrame.x,xbins);
[positions_per_imageFrame.yd] = discretize(positions_per_imageFrame.y,ybins);

place_cell_map = TNC_CreateRBColormap(1024,'exag');
pc_kern = TNC_CreateGaussian(100,2,200,1);

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

    figure(102); clf;
    for zz = 1:numel(inds)
        subplot(13,12,zz);
        tmp = mean_resp(:,:,zz);
        tmp(isnan(tmp)) = 0;
        imagesc(xbins,ybins,conv2(tmp,pc_kern'*pc_kern,'same'),[0 7]); colormap(place_cell_map); axis off; title(num2str(max(conv2(tmp,pc_kern'*pc_kern,'same'),[],'all')))
        drawnow;
    end

    %% Try to find sequences of activity

    ord_map = TNC_CreateRBColormap(31,'cpob');

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
    clust_inds = perm(90:120)

    figure(2); clf;
    for vv=clust_inds

        plot(valid_frames,trace_mat.dffxPos(vv,:)+(10*find(vv==clust_inds)),'color',ord_map(find(vv==clust_inds),:)); hold on; xlim([0 2e4]);

    end

    figure(3); clf;
    plot(valid_frames,positions_per_imageFrame.x,'k','Color',[0 0 0 1]); xlim([0 2e4]);
    hold on;
    plot(valid_frames,positions_per_imageFrame.y,'k','Color',[1 0 0 1]); xlim([0 2e4]);
    plot(valid_frames,5*sum(trace_mat.dffxPos(perm,:),1),'k','Color',[0 0.67 1 1]); axis tight; xlim([0 2e4]); 
    