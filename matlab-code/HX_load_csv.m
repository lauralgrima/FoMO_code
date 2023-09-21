function [hexa_data] = HX_load_csv(filename, verbose, photo)

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

if photo
    TP = readtable([filename(1:end-7) 'photo.csv']);
    hexa_data.photo.dFF  = TP.signal;
    hexa_data.photo.sess = TP.session_number;
end
if verbose
    hexa_data
end