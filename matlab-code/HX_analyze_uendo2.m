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

% can be aligned by comparing to tracker.Ext_ClockCount which increments by
% 1 each time a GPIO pulse is received

%% New loading the compiled data

traces          = readtable('traces_tracking.csv');
behav           = readtable('beh.csv');

