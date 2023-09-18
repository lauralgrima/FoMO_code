% HX_analyze_microendo

%% Data Loading / Construction of data matrix
clear trace_mat

filename_traces = '6PM4_2022-10-14_celltraces.csv';
filename_traces_props = '6PM4_2022-10-14_celltraces-props.csv';
filename_spks = '6PM4_2022-10-14_spikes.csv';
filename_spks_props = '6PM4_2022-10-14_spikes-props.csv';

exag = TNC_CreateRBColormap(1000,'exag');
catmap = TNC_CreateRBColormap(1000,'cat2');

traces = readtable(filename_traces);
traces_props = readtable(filename_traces_props);
spks = readtable(filename_spks);
spks_props = readtable(filename_spks_props);

trace_mat.dff = [];

for kk=1:numel(traces.Properties.VariableNames)
    
    if numel(strfind(traces.Properties.VariableNames{kk},'C'))>0
        tmp = traces.(traces.Properties.VariableNames{kk});
        trace_mat.dff = [trace_mat.dff ; tmp'];
    elseif numel(strfind(traces.Properties.VariableNames{kk},'Var1'))
        trace_mat.t = traces.(traces.Properties.VariableNames{kk});
    end
    
end


for kk=1:numel(traces.Properties.VariableNames)
    
    if numel(strfind(traces.Properties.VariableNames{kk},'C'))>0
        
        event_inds = contains(spks.CellName,traces.Properties.VariableNames{kk});
        spks.Time_s_(event_inds);
        spk_mat.cell(str2num(traces.Properties.VariableNames{kk}(2:end))+1).t = spks.Time_s_(event_inds);
        spk_mat.cell(str2num(traces.Properties.VariableNames{kk}(2:end))+1).v = spks.Value(event_inds);
        
    end
end

trace_mat.x = traces_props.CentroidX;
trace_mat.y = traces_props.CentroidY;


mappedX = tsne(trace_mat.dff);

[idx,C,sumdist] = kmeans(mappedX,5,'Distance','cityblock', 'Display','final','Replicates',10);
[idxs,inds]=sort(idx);

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


spk_mat.log = zeros(size(trace_mat.dff));
spk_mat.cnt = zeros(size(trace_mat.dff));
figure(10); clf;
subplot(5,1,1:4);

for kk=1:numel(spk_mat.cell)
    
    for mm=1:numel(spk_mat.cell(kk).t)
        spk_mat.cell(kk).c(mm) = find(trace_mat.t<=spk_mat.cell(kk).t(mm),1,'last');
    end
    
    spk_mat.log(kk,spk_mat.cell(kk).c) = 1;
    spk_mat.cnt(kk,spk_mat.cell(kk).c) = spk_mat.cell(kk).v;
    
    scatter(spk_mat.cell(kk).t,find(inds==kk).*ones(1,numel(spk_mat.cell(kk).t)),1,catmap(idx(kk),:),'|','LineWidth',2); hold on; colormap(catmap);
    
end

subplot(5,1,5);
spe_kern = [0 ones(1,5) 0]./5;
spe_kern = TNC_CreateGaussian(500,5,1000,1); spe_kern = spe_kern./max(spe_kern);
plot([1:size(spk_mat.log,2)].*0.05,conv(sum(spk_mat.log,1),spe_kern,'same'));
box off;


%% Load behavioral tracking data and rewards etc and sort out alignment
debug = 0;
filename_track = '../raw/6PM4_hex_foraging_V2_2022-10-14-102642.csv';
tracker = readtable(filename_track);
filename_behav = '6PM4_behaviour.csv';
behaviour = readtable(filename_behav);

behav.pos.y = tracker.YLocation .* 0.3;
fill = find( behav.pos.y == 0 );
cont = find( behav.pos.y ~= 0 );

behav.pos.y(fill) = interp1(cont,behav.pos.y(cont),fill);

tmp_img_x = tracker.XLocation.* 0.3;
tmp_img_x(fill) = interp1(cont,tmp_img_x(cont),fill);

behav.pos.x = tracker.StagePostion.* 0.0069 + tmp_img_x;

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
end

figure(101); clf;
plot(sgolayfilt( behav.pos.x , 3 , 101) , sgolayfilt( behav.pos.y , 3, 101) ); axis([0 2500 0 500]);
axis equal; grid on; box off;


% rewards

rew_inds = find(behaviour.rewarded==1);
behaviour.event_time(rew_inds);
figure(102); scatter(behaviour.event_time(rew_inds),ones(1,numel(rew_inds)),100,'|');

behav.rewards_micro_i = behaviour.micro_i(rew_inds);

figure(10); 
subplot(5,1,5);
hold on;
scatter(behav.rewards_micro_i,ones(1,numel(behav.rewards_micro_i)),50,'r','|');

figure(11);
[sink] = TNC_ExtTrigWins(conv(sum(spk_mat.log,1),spe_kern,'same'),behav.rewards_micro_i.*20,[60 180]);
shadedErrorBar(0.05*sink.range,sink.avg,sink.err);
ylabel('Synchronous activity'); xlabel('Time from reward (s)'); box off; grid on;

%%
filename_gpio = '../raw/6PM4_2022-10-14_GPIO.csv';
gpio = readtable(filename_gpio);
% clear behav;

tf = contains(gpio.ChannelName,'GPIO-1');
behav.val = gpio.Value(tf==1)';
times = gpio.Time_s_(tf==1)';

pulses = find(behav.val>1e4 & [0 diff(behav.val)]>0);

figure(8); clf; plot(behav.val);
hold on;
plot(pulses,ones(1,numel(pulses)).*1e4,'r*');
start_pulses = pulses(find([251 diff(pulses)>250]));
plot(start_pulses,ones(1,numel(start_pulses)).*1e4,'go');

behav.gpio.times = times(start_pulses);

% can be aligned by comparing to tracker.Ext_ClockCount which increments by
% 1 each time a GPIO pulse is received


