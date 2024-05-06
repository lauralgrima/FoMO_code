%% TNC_S_6PortForager

% % DEPENDS UPON https://github.com/dudmanj/TONIC

% Notes on metadata files:
% The columns in the csv are:
% 1. Frame no.
% 2. Rsync pulse count (cumulative)
% 3. x position
% 4. y position
% 5. Stage position
% X and Y are centroids of the tracking box

%% LOAD csv data

fname = 'Metadata_210302130619_6P9.csv';
% cent_kern = TNC_CreateGaussian(50,2,100,1);
cent_kern = [0 1 1 1 0]/3;
% close all;
session = dlmread(fname);
win = 5;
start = 500;

%seems to be a y offset of something like ~50 pixels
y_off=80;

% img_x = sgolayfilt( (session(start:end,3)-median(session(start:end,3))) .* 0.3 , 3 , 21); % scaled into mm
pxl_x = session(start:end,3);
    pxl_x(pxl_x>320 & pxl_x<1280) = 320;
    pxl_x(pxl_x>=1280) = pxl_x(pxl_x>=1280)-1280+320;
%     pxl_x(pxl_x<=320) = 320 + (pxl_x(pxl_x<=320)-320);
pxl_y = session(start:end,4)+y_off; % scaled into mm
    pxl_y(pxl_y>320 & pxl_y<1280) = 320;
    pxl_y(pxl_y>=1280) = pxl_y(pxl_y>=1280)-1280+320;
%     pxl_y(pxl_y<=320) = pxl_y(pxl_y<=320)+min(pxl_y);

% Original version trying to just smooth/filter
img_x = sgolayfilt( (session(start:end,3)-median(session(start:end,3))) .* 0.3 , 3, 151);
img_y = sgolayfilt( (session(start:end,4)-min(session(start:end,3))) .* 0.3 , 3, 151); % scaled into mm
stg_x = session(start:end,5) .* 0.0069; % scaled into mm
frms = session(start:end,1)';

% pos.x = img_x+stg_x;
% % pos.x = stg_x;
% pos.y = img_y;
% 
% figure(1); clf; plot(pos.xf,pos.yf,'color',[0.5 0.5 0.5 0.1]); %hold on; plot(pos.x(tmp2),pos.y(tmp2),'c*'); axis off; axis tight;
% 
% figure(2); clf; plot(sgolayfilt(conv(pos.v,[0 ones(1,25)/25 0],'same'),3,11),'color',[0.5 0.5 0.5]); hold on; plot(conv(sqrt([0 ; diff(stg_x).^2]),[0 ones(1,25)/25 0],'same'),'color',[0.5 0 0]); plot(conv(sqrt([0 ; diff(img_y).^2]),[0 ones(1,25)/25 0],'same'),'color',[0 0.33 0.5]); 


%% Convert stage movements to smooth position (interpolate over encoder reads)

stg_mv = find(abs([0 diff(stg_x')])>0);
[x_stg_align] = TNC_ExtTrigWins(img_x,stg_mv,[10 10]);

figure(9); clf; plot([0 diff(stg_x')]); hold on; plot(stg_x'); plot(stg_mv,stg_x(stg_mv),'o');

stg_mvL = abs([0 diff(stg_x')])>0;
stg_xi = stg_x';
stg_xi(stg_mv) = stg_x(stg_mv)';
stg_xi(~stg_mvL) = interp1(frms(stg_mvL),stg_x(stg_mvL)',frms(~stg_mvL));

figure(10); clf; plot(stg_x'); hold on; plot(stg_xi);
pos.x = stg_xi+img_x';
pos.y = img_y';

% figure(10); clf; plot(x_stg_align.range,x_stg_align.avg);


%% Convert to position and speed measurements
tmp = find(img_x==0 & img_y==0);
not_tmp = find(img_x~=0 & img_y~=0);
pos.x(tmp) = interp1(not_tmp,pos.x(not_tmp),tmp);
pos.y(tmp) = interp1(not_tmp,pos.y(not_tmp),tmp);

pos.v = sqrt( [0 , diff(pos.x)].^2 + [0 , diff(pos.y)].^2 );
tmp2 = find(pos.v>90);
not_tmp2 = find(pos.v<=0);
if numel(not_tmp2) > 1
    pos.v(tmp2) = interp1(not_tmp2,pos.v(not_tmp2),tmp2);
end

pos.xf = sgolayfilt( pos.x , 3 , 31 );
pos.yf = sgolayfilt( pos.y , 3 , 31 );


%% Extract run events

pos.v_s = sgolayfilt(conv(pos.v,[0 ones(1,25)/25 0],'same'),3,11);
[events] = TNC_ExtractPeaks(pos.v_s,1.75,150,1);

% Need to think about merging segments here

% plot out run trajectories
figure(5); clf; 
% plot(pos.xf,pos.yf,'color',[0.5 0.5 0.5 0.1]); hold on;

clear runs; cnt = 1;
dist_thresh = 200;
dist = zeros(1,events.num);

for pp=1:events.num
%     clf; plot(pos.xf,pos.yf,'color',[0.5 0.5 0.5 0.1]);  axis tight; axis off;

    % compute distance of run
    dist(pp) = abs(pos.xf(events.stops(pp))-pos.xf(events.starts(pp))) + abs(pos.yf(events.stops(pp))-pos.yf(events.starts(pp)));
    
    % threshold distance and update cleaned runs
    if dist(pp)>dist_thresh
        
        runs.starts(cnt) = events.starts(pp);
        runs.stops(cnt) = events.stops(pp);
        runs.dist(cnt) = dist(pp);
        
        % if passes distance threshold then also plot
        if rand(1)<0.05
            plot(pos.xf(events.starts(pp):events.stops(pp)),pos.yf(events.starts(pp):events.stops(pp)),'color',[0 0.33 0.5 0.5],'linewidth',3); 
        else
            plot(pos.xf(events.starts(pp):events.stops(pp)),pos.yf(events.starts(pp):events.stops(pp)),'color',[0 0.33 0.5 0.05],'linewidth',1); 
        end
        hold on; 
        plot(pos.xf(events.starts(pp)),pos.yf(events.starts(pp)),'*','color',[0 1 0],'linewidth',2);
        plot(pos.xf(events.stops(pp)),pos.yf(events.stops(pp)),'*','color',[1 0 0],'linewidth',2);
        if cnt==1
            axis([-250 2000 -15 435]);
            grid on; box off;
            pbaspect([5 1 1]);
            xlabel('mm');
            ylabel('mm');
            set(gca,'YDir','reverse');
        end
%         drawnow;
        cnt = cnt+1;
    end
end

%% Examine runs by extracting aligned video frames

[~,eg_traj] = max(runs.dist);
eg_traj = 1;
writeMov = 1;

target_v = VideoReader('Video_210302130619.avi');

if writeMov
    v = VideoWriter(['ExampleTrajectory_' num2str(eg_traj) '.mp4'],'MPEG-4');
    v.Quality = 95;
    v.FrameRate = 150;
    open(v);
end

for pp=eg_traj %1:numel(runs.dist)
        
    start_frm = runs.starts(pp);
    stop_frm = runs.stops(pp);
    
    start_frm = 11000;
    stop_frm = 31000;

    frm_vec = start_frm:stop_frm;
    
    figure(8); clf;
    for kk=frm_vec
        subplot(1,10,1:8);
        hold off;
        plot(pos.xf(start_frm:stop_frm),pos.yf(start_frm:stop_frm),'color',[0 0.33 0.5 0.5],'linewidth',1); hold on;
        plot(pos.xf(kk),pos.yf(kk),'.','MarkerSize',50,'color',[0.5 0 0 0.5],'linewidth',5);
        axis([-250 2000 -15 435]);
        set(gca,'YDir','reverse');
        grid on; box off;
        pbaspect([5 1 1]);
        xlabel('mm');
        ylabel('mm');

        subplot(1,10,9:10); hold off;
        frame = rgb2gray(read(target_v,frms(kk)));
        imshow(imadjust(frame)); axis off;
        pbaspect([1 1 1]);
        hold on;
        plot(pxl_x(kk),pxl_y(kk),'r+','linewidth',5);
        drawnow;
        
        if writeMov
            frameW = getframe(gcf);
            writeVideo(v,frameW);
        end
        
    end
end

if writeMov
    close(v);
end


%% Lana DBS data analysis

clear pos;

true_pos_ts = find(~isnan(dbs5.x_loc));

pos.x = dbs5.x_loc(true_pos_ts);
pos.y = dbs5.y_loc(true_pos_ts);

true_pos_ph = find(~isnan(dbs5.photo_red));

photo.r = dbs5.photo_red(true_pos_ph);
photo.g = dbs5.photo_green(true_pos_ph);
photo.c = dbs5.photo_corr(true_pos_ph);

% correct tracking errors
trk_lost = find(pos.x==0);
trk_true = find(pos.x~=0);
pos.xi = pos.x;
pos.yi = pos.y;

pos.xi(trk_lost) = interp1(trk_true,pos.xi(trk_true),trk_lost);
pos.yi(trk_lost) = interp1(trk_true,pos.yi(trk_true),trk_lost);

pos.xs = sgolayfilt(pos.xi,3,21);
pos.ys = sgolayfilt(pos.yi,3,21);
pos.spd = sgolayfilt( sqrt( [0 ; diff(pos.xs)].^2 + [0 ; diff(pos.ys)].^2 ) , 3, 21);
pos.spd(pos.spd<=0) = 0;
pos.spdL = conv( pos.spd ,[0 ones(1,500)/500 0],'same' );


dbs_on = find(~isnan(dbs5.rsync));
for gg=1:numel(dbs_on)
    dbs_on_pos(gg) = find(true_pos_ts>dbs_on(gg),1);
    dbs_off_pos(gg) = find(true_pos_ts>dbs_on(gg)+60e3,1);
end

[dbs_spd] = TNC_ExtTrigWins(pos.spdL,dbs_on_pos,[3000,15000]);
[dbs_xs] = TNC_ExtTrigWins(pos.xs,dbs_on_pos,[3000,15000]);
[dbs_ys] = TNC_ExtTrigWins(pos.ys,dbs_on_pos,[3000,15000]);

shuffs = 1000;
shuff_data = zeros(shuffs,numel(dbs_spd.range));
non_dbs_times = [ ];
for qq=1:numel(dbs_off_pos)-1
    non_dbs_times = [  non_dbs_times dbs_off_pos(qq)+15e3:dbs_on_pos(qq+1)-15e3 ];
end

for qq=1:1000
    [shf_spd] = TNC_ExtTrigWins(pos.spdL,non_dbs_times(randperm(numel(non_dbs_times),5)),[3000,15000]);
    shf_data(qq,:) = shf_spd.avg - mean(shf_spd.avg(1:3000));
end

figure(100); clf;
subplot(121);
plot(pos.xs,pos.ys,'color',[0.8 0.8 0.8]); hold on;
for yy = 1:numel(dbs_on_pos)
plot(pos.xs(dbs_on_pos(yy):dbs_off_pos(yy)),pos.ys(dbs_on_pos(yy):dbs_off_pos(yy)));
end

subplot(122);
plot(pos.spdL,'color',[0.8 0.8 0.8]); hold on;
for yy = 1:numel(dbs_on_pos)
plot(dbs_on_pos(yy):dbs_off_pos(yy),pos.spdL(dbs_on_pos(yy):dbs_off_pos(yy)));
end
% yyaxis right;
% plot(photo.r); hold on;
% plot(photo.g);
% plot(photo.c);

figure(101); clf;
shadedErrorBar(dbs_spd.range,mean(shf_data),std(shf_data),{'color',[0 0 0 0.1]});
hold on;
plot(dbs_spd.range,dbs_spd.avg-mean(dbs_spd.avg(1:3000)),'linewidth',2,'color',[0.7 0 0]);
box off;
ylabel('Baseline subtracted speed'); xlabel('Time (ms)');

%% Lana DBS data analysis

clear pos dbs_*;

true_pos_ts = find(~isnan(dbs7.x_loc));

pos.x = dbs7.x_loc(true_pos_ts);
pos.y = dbs7.y_loc(true_pos_ts);


% correct tracking errors
trk_lost = find(pos.x==0);
trk_true = find(pos.x~=0);
pos.xi = pos.x;
pos.yi = pos.y;

pos.xi(trk_lost) = interp1(trk_true,pos.xi(trk_true),trk_lost);
pos.yi(trk_lost) = interp1(trk_true,pos.yi(trk_true),trk_lost);

pos.xs = sgolayfilt(pos.xi,3,21);
pos.ys = sgolayfilt(pos.yi,3,21);
pos.spd = sgolayfilt( sqrt( [0 ; diff(pos.xs)].^2 + [0 ; diff(pos.ys)].^2 ) , 3, 21);
pos.spd(pos.spd<=0) = 0;
pos.spdL = conv( pos.spd ,[0 ones(1,1000)/1000 0],'same' );


dbs_on = find(~isnan(dbs7.rsync));
for gg=1:numel(dbs_on)
    dbs_on_pos(gg) = find(true_pos_ts>dbs_on(gg),1);
    dbs_off_pos(gg) = find(true_pos_ts>dbs_on(gg)+60e3,1);
end

[dbs_spd] = TNC_ExtTrigWins(pos.spdL,dbs_on_pos,[3000,15000]);

shuffs = 1000;
shuff_data = zeros(shuffs,numel(dbs_spd.range));
non_dbs_times = [ ];
for qq=1:numel(dbs_off_pos)-1
    non_dbs_times = [  non_dbs_times dbs_off_pos(qq)+15e3:dbs_on_pos(qq+1)-15e3 ];
end

for qq=1:1000
    [shf_spd] = TNC_ExtTrigWins(pos.spdL,non_dbs_times(randperm(numel(non_dbs_times),5)),[3000,15000]);
    shf_data(qq,:) = shf_spd.avg - mean(shf_spd.avg(1:3000));
end

figure(110); clf;
subplot(121);
plot(pos.xs,pos.ys);

subplot(122);
plot(pos.spd);

figure(111); clf;
shadedErrorBar(dbs_spd.range,mean(shf_data),std(shf_data),{'color',[0 0 0 0.1]});
hold on;
plot(dbs_spd.range,dbs_spd.avg-mean(dbs_spd.avg(1:3000)),'linewidth',2,'color',[0.7 0 0]);
box off;
ylabel('Baseline subtracted speed'); xlabel('Time (ms)');

