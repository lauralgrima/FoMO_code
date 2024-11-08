function [pos] = TNC_BigBoiConverter(session)

cent_kern = [0 1 1 1 0]/3;
win = 5;
start = find(session(:,3)>0 & session(:,6)>0,1);
frms = session(start:end,1)';

% Convert stage movements to smooth position (interpolate over encoder reads)

pxl_x = sgolayfilt( session(start:end,3).* 0.3, 3, 31);
pxl_y = sgolayfilt( session(start:end,6).* 0.3, 3, 31);
stg_x = session(start:end,7).* 0.0069;

pos.x = stg_x'+pxl_x';
pos.y = pxl_y';

figure(10); clf; 
plot(pos.x,pos.y,'color',[0 0 0 0.1]);

% Convert to position and speed measurements
tmp = find(pxl_x==0 & pxl_y==0);
not_tmp = find(pxl_x~=0 & pxl_y~=0);
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
pos.frms = frms;

figure(10); clf; plot(pos.xf,pos.yf,'color',[0 0 0 0.1]); box off; axis off;

%% Deprecated code from trying to get rid of small stage move x artefact
% appears unnecessary in newer files.

% stg_mv = find(abs([0 diff(stg_x')])>0);
% stg_mv = stg_mv(stg_mv>10 & stg_mv<numel(pxl_x)-10);
% [x_stg_align] = TNC_ExtTrigWins(pxl_x,stg_mv,[10 10]);
% stg_mvL = abs([0 diff(stg_x')])>0;
% stg_xi = stg_x';
% stg_xi(stg_mv) = stg_x(stg_mv)';
% stg_xi(~stg_mvL) = interp1(frms(stg_mvL),stg_x(stg_mvL)',frms(~stg_mvL));