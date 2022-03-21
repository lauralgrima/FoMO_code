%% Notes on the hexaport hackathon

% function 
% Pass in data file and policy choice
% Format the simulation to produce same size output
% Plot and compare model to data

max_tsteps = 24000 * 10; % 24000 seconds * 10 Hz = 10 rewards at worst baited port

for t=1:max_tsteps
    
   % should we check any port at this time point
   
   % discount reward probabilities by distance to port
   
   % which port should we check? softmax{ Pr(R|port) }

   % Update our 'policy' Pr(R|port)
   switch policy
       
       case ''
      
       otherwise % random policy
           
   end
    
    
end

%% Exploring the licking data a bit

lk_kern = TNC_CreateGaussian(500,200,1000,1);
lk_kern = lk_kern./ max(lk_kern);
% lk_kern(1:500)=0;
% lk_kern = [0 ones(1,1000) 0];

% make a matrix of dimension 6 ports x total time
lk_6p6_raw = csvread('lick_df.csv',1,1);
hexa.all_raw = zeros(6,lk_6p6_raw(end,2));
prior_port_ids = [lk_6p6_raw(1,3) ; lk_6p6_raw(1:end-1,3)]';
for qq=1:6
    hexa.all_raw(qq,lk_6p6_raw(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1,1)) = 1;
    % find the prior port
    hexa.prior_port(qq,lk_6p6_raw(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1,1)) = prior_port_ids(lk_6p6_raw(:,3)==qq & lk_6p6_raw(:,2)==1);
end

figure(11); clf;
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
end

eg_port = 2;
figure(12); clf;
subplot(2,4,1:2);
for qq=1:6
    plot(hexa.all_raw(qq,:)+qq); hold on;
end
box off;

lk_intervals = [0 diff(find(hexa.all_raw(eg_port,:)==1))];
lk_times = find(hexa.all_raw(eg_port,:)==1);

subplot(2,4,5:6);
plot(hexa.all_dat(eg_port,:)); hold on;
plot(conv(hexa.all_raw(eg_port,:),lk_kern,'same'));
plot(lk_times(lk_intervals>1000),ones(1,sum(lk_intervals>1000)),'*');
box off;

subplot(2,4,[3 7]);
hist(log10(diff(find(hexa.all_raw(eg_port,:)==1))),250);


% add a piece of logic to check for switched port
port_ids = hexa.all_raw .* [1:6]';

figure(13); clf;
for qq=1:6
    lk_intervals = [0 diff(find(hexa.all_raw(qq,:)==1))];
    port_changes = [hexa.prior_port(qq,find(hexa.all_raw(qq,:)==1))];
    visit_cnts(qq,1) = sum(lk_intervals>1000);
    visit_cnts(qq,2) = sum(port_changes~=qq);
    subplot(1,6,qq);
    swarmchart(port_changes,log10(lk_intervals),'.'); hold on;
    plot([0 7],[3 3],'k-');
    axis([0 7 0 7]);
    xlabel('Previous port'); ylabel('Lick interval (log ms)'); title(['port: ' num2str(qq)]);
end
figure(12);
subplot(2,4,[4 8]);
plot(1:6,sum(hexa.all_raw,2),'o-');
yyaxis right
plot(1:6,visit_cnts(:,1),'o-');
plot(1:6,visit_cnts(:,2),'o-');

