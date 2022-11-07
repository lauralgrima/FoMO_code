cd('/Users/grimal/Dropbox (HHMI)/hexaport/hdf5');

metadata  = h5readatt('6PG5_NAc_conc_behav.h5','/mydata','metadata');
animal_no = metadata((17:20));
region    = metadata((30:35));
task      = metadata((60:63));
% can also include check threshold, number of sessions

% show what hdf5 contains
h5disp('6PG5_NAc_conc_behav.h5')

columns = h5read('6PG5_NAc_conc_behav.h5','/mydata/axis0');
index = h5read('6PG5_NAc_conc_behav.h5','/mydata/axis1');
dates = h5read('6PG5_NAc_conc_behav.h5','/mydata/block0_values');
stuff = h5read('6PG5_NAc_conc_behav.h5','/mydata/block2_values');
session_n = stuff(1,:); % 1 to 5
unique_vis = stuff(2,:); % unique port visit

everything_else = h5read('6PG5_NAc_conc_behav.h5','/mydata/block1_values'); % rest of the data
event_time      = everything_else(1,:); % in minutes
on_off          = everything_else(2,:); % only lick on included (always 1 or NaN)
port_n          = everything_else(3,:);
port_rank       = everything_else(4,:);
prior_port      = everything_else(5,:);
rewarded        = everything_else(6,:);
RA_1            = everything_else(7,:); % RA is one when reward becomes available at the port 
RA_2            = everything_else(8,:);
RA_3            = everything_else(9,:);
RA_4            = everything_else(10,:);
RA_5            = everything_else(11,:);
RA_6            = everything_else(12,:);
photo_i         = everything_else(13,:); % for aligning with photometry
video_i         = everything_else(14,:); % for aligning with video frames



