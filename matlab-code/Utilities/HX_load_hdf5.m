function [hexa_data] = HX_load_hdf5(path , filename, verbose)

metadata  = h5readatt([path '/' filename],'/mydata','metadata');
hexa_data.animal_no = metadata([14:17]+strfind(metadata,'animal_no'));
hexa_data.region    = metadata([11:13]+strfind(metadata,'region'));
hexa_data.task      = metadata([9:12]+strfind(metadata,'task'));
% can also include check threshold, number of sessions

if verbose
    % show what hdf5 contains
    h5disp([path '/' filename])
end
hexa_data.filename=filename;

hexa_data.columns = h5read([path '/' filename],'/mydata/axis0');
hexa_data.index = h5read([path '/' filename],'/mydata/axis1')';
hexa_data.dates = h5read([path '/' filename],'/mydata/block0_values');
stuff = h5read([path '/' filename],'/mydata/block2_values');
hexa_data.session_n = stuff(1,:); % 1 to 5
hexa_data.unique_vis = stuff(2,:); % unique port visit

everything_else = h5read([path '/' filename],'/mydata/block1_values'); % rest of the data
hexa_data.event_time      = everything_else(1,:); % in minutes
hexa_data.on_off          = everything_else(2,:); % only lick on included (always 1 or NaN)
hexa_data.port_n          = everything_else(3,:);
hexa_data.port_rank       = everything_else(4,:);
hexa_data.prior_port      = everything_else(5,:);
hexa_data.rewarded        = everything_else(6,:);
hexa_data.RA_1            = everything_else(7,:); % RA is one when reward becomes available at the port 
hexa_data.RA_2            = everything_else(8,:);
hexa_data.RA_3            = everything_else(9,:);
hexa_data.RA_4            = everything_else(10,:);
hexa_data.RA_5            = everything_else(11,:);
hexa_data.RA_6            = everything_else(12,:);
hexa_data.photo_i         = everything_else(13,:); % for aligning with photometry
hexa_data.video_i         = everything_else(14,:); % for aligning with video frames



