% HX_analyze_microendo
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

rew_inds = find(behaviour.unique_visit==1 & behaviour.rewarded==1);

behaviour.event_time(rew_inds);

figure(102); scatter(behaviour.event_time(rew_inds),ones(1,numel(rew_inds)),100,'|');

behav.rewards_micro_i = behaviour.micro_i(rew_inds);
behav.rewards_port_id = behaviour.port(rew_inds);

visit_inds          = find(behaviour.unique_visit==1 & behaviour.rewarded==0);

behav.visit_port_id = behaviour.port(visit_inds);
behav.visit_micro_i = behaviour.micro_i(visit_inds);


figure(10); clf;
scatter(behav.rewards_micro_i,behav.rewards_port_id+30,50,behav.rewards_port_id,'|','linewidth',3);
hold on;
scatter(behav.visit_micro_i,behav.visit_port_id+30,50,behav.visit_port_id,'|','linewidth',1);

figure(11); clf;
[sink] = TNC_ExtTrigWins(conv(sum(spk_mat.log,1),spe_kern,'same'),behav.rewards_micro_i.*20,[60 120]);
shadedErrorBar(0.05*sink.range,sink.avg,sink.err,{'color',[0 0.33 0.55]});
ylabel('Synchronous activity'); box off; grid on;
hold on;
[sink] = TNC_ExtTrigWins(conv(sum(spk_mat.log,1),spe_kern,'same'),behav.visit_micro_i.*20,[60 120]);
shadedErrorBar(0.05*sink.range,sink.avg,sink.err);
ylabel('Synchronous activity'); xlabel('Time from visit (s)'); box off; grid on;


[pv_rew] = TNC_ExtTrigWins3d(spk_mat.cnt,behav.rewards_micro_i.*20,[60 120]);
for zz = 1:size(pv_rew.win_avg,1)
    [~,lag] = max(conv(pv_rew.win_avg(zz,:),spe_kern,'same'));
    for_sort_cells(zz) = lag;
end

[~,c_inds] = sort(for_sort_cells,'ascend');
figure(19); clf; imagesc(pv_rew.win_avg(c_inds,:),[0 1]); colormap(exag);

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

figure(10); hold on;

for kk=1:numel(spk_mat.cell)
    
    for mm=1:numel(spk_mat.cell(kk).t)
        spk_mat.cell(kk).c(mm) = find(trace_mat.t<=spk_mat.cell(kk).t(mm),1,'last');
    end
    
    spk_mat.log(kk,spk_mat.cell(kk).c) = 1;
    spk_mat.cnt(kk,spk_mat.cell(kk).c) = spk_mat.cell(kk).v;
    
    scatter(spk_mat.cell(kk).t,50+(find(inds==kk).*ones(1,numel(spk_mat.cell(kk).t))),1,catmap(idx(kk)+2,:),'|','LineWidth',2); hold on; colormap(catmap);
    
%     scatter(spk_mat.cell(kk).t,50+(find(c_inds==kk).*ones(1,numel(spk_mat.cell(kk).t))),1,catmap(idx(kk)+2,:),'|','LineWidth',2); hold on; colormap(catmap);
    
end

hold on;
spe_kern = [0 ones(1,5) 0]./5;
spe_kern = TNC_CreateGaussian(500,5,1000,1); spe_kern = spe_kern./max(spe_kern);
plot([1:size(spk_mat.log,2)].*0.05,conv(sum(spk_mat.log,1),spe_kern,'same'),'color',[0 0 0 0.5]);
box off;

%% Get x and y pos at points along the spk_mat
valid_samps = find(~isnan(behaviour.micro_i));
P = polyfit(behaviour.micro_i(valid_samps),behaviour.vid_i(valid_samps),1);

if debug
    figure(200); scatter(behaviour.micro_i(valid_samps),behaviour.vid_i(valid_samps));
end

tracker_pnts = round(trace_mat.t * P(1))+1;

% pretty sure there is a minor offset issue here - matches the estimate
% from linear fit, so zero padding to get offset assuming stable sample
% rate
spk_mat.pos.x = [zeros(1,531) decimate(behav.pos.x,10)'];
spk_mat.pos.y = [zeros(1,531) decimate(behav.pos.y,10)'];

%%
figure(201); clf;
example_cells = randperm(144,60);
cnt = 1;



for qq=example_cells

    subplot(10,6,cnt);
    plot(spk_mat.pos.x,spk_mat.pos.y,'color',[0 0 0 0.1]); hold on;
    spk_times = find(spk_mat.log(qq,:)~=0);
    scatter(spk_mat.pos.x(spk_times(spk_times>531)),spk_mat.pos.y(spk_times(spk_times>531)),10,'r','filled');
    axis off; axis equal; axis off;
    cnt=cnt+1;
end

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


