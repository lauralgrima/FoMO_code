function [hexa_data] = HX_load_csv(filename, verbose, photo, photo_filename)

T = readtable(filename);
hexa_data.filename=filename;
tmp = strfind(filename,'_');

hexa_data.animal_no = filename(1:tmp(1)-1);
hexa_data.region    = filename(tmp(1)+[1:3]);
hexa_data.task      = filename(tmp(2)+1:tmp(3)-1);
hexa_data.dates     = T.date;
hexa_data.session_n = T.session_number; % 1 to 5
hexa_data.unique_vis      = T.unique_visit; % unique port visit
hexa_data.event_time      = T.event_time; % in minutes
hexa_data.on_off          = T.on_off; % only lick on included (always 1 or NaN)
hexa_data.port_n          = T.port;
hexa_data.port_rank       = T.port_rank;
hexa_data.rewarded        = T.rewarded;
hexa_data.RA_1            = T.RA_1; % RA is an integer reflecting when reward becomes available; value is total available to that point
hexa_data.RA_2            = T.RA_2;
hexa_data.RA_3            = T.RA_3;
hexa_data.RA_4            = T.RA_4;
hexa_data.RA_5            = T.RA_5;
hexa_data.RA_6            = T.RA_6;
hexa_data.photo_i         = T.photo_i; % for aligning with photometry
hexa_data.video_i         = T.video_i; % for aligning with video frames

hexa_data.photo_i_con     = hexa_data.photo_i(hexa_data.session_n==1);
hexa_data.event_time_con  = hexa_data.event_time(hexa_data.session_n==1);
shift = 0; shift_e = 0;

if photo

    TP = readtable(photo_filename);
    hexa_data.photo.dFF  = TP.signal;
    hexa_data.photo.sess = TP.session_number;
end

for zz=2:max(unique(hexa_data.session_n))
    shift                   = find(hexa_data.photo.sess==zz-1,1,'last');
    hexa_data.photo_i_con   = [hexa_data.photo_i_con ; hexa_data.photo_i(hexa_data.session_n==zz)+shift];
    shift_e                 = shift_e+max(hexa_data.event_time(hexa_data.session_n==zz-1));
    hexa_data.event_time_con= [hexa_data.event_time_con ; hexa_data.event_time(hexa_data.session_n==zz)+shift_e];
end


if verbose
    hexa_data
end