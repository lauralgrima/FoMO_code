function [hexa_data_an] = HX_analyze_session(hexa_data,session)

all_rew_inds = find(hexa_data.rewarded==1 & hexa_data.session_n==session);
all_rew_ports = hexa_data.port_n(all_rew_inds);

all_vis_inds = find(hexa_data.unique_vis==1 & hexa_data.session_n==session);
all_vis_ports = hexa_data.port_n(all_vis_inds);

figure(57); clf;
cat_map = TNC_CreateRBColormap(8,'mapb');

hexa_data_an.visits = zeros(6,max(find(hexa_data.session_n==session)));
hexa_data_an.rewards = zeros(6,max(find(hexa_data.session_n==session)));

hexa_data_an.vi.avg = zeros(1,6);
hexa_data_an.vi.std = zeros(1,6);
hexa_data_an.rw.avg = zeros(1,6);
hexa_data_an.rw.std = zeros(1,6);

for qq=unique(all_rew_ports)

    subplot(131)
    hexa_data_an.port_rew(qq).ts = all_rew_inds(all_rew_ports==qq);
    hexa_data_an.port_rew(qq).iti = diff(hexa_data_an.port_rew(qq).ts);
    disp(['Port ' num2str(qq) ' mean ri: ' num2str(mean(hexa_data_an.port_rew(qq).iti)) ' +/- ' num2str(std(hexa_data_an.port_rew(qq).iti))]);
    hexa_data_an.port_rew(qq).histX = 0.1:0.1:4;
    hexa_data_an.port_rew(qq).histY = hist(log10(hexa_data_an.port_rew(qq).iti),hexa_data_an.port_rew(qq).histX);
    plot(hexa_data_an.port_rew(qq).histX , cumsum(hexa_data_an.port_rew(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','2','3','4','5','6'},'Location','northwest');
        xlabel('Log10 Interreward Interval'); ylabel('Count'); box off;
    end

    hexa_data_an.rewards(qq,all_rew_inds(all_rew_ports==qq)) = 1;
    hexa_data_an.rw.avg(qq) = mean(hexa_data_an.port_rew(qq).iti);
    hexa_data_an.rw.std(qq) = std(hexa_data_an.port_rew(qq).iti);

    subplot(132)
    hexa_data_an.port_vis(qq).ts = all_vis_inds(all_vis_ports==qq);
    hexa_data_an.port_vis(qq).iti = diff(hexa_data_an.port_vis(qq).ts);
    disp(['Port ' num2str(qq) ' mean vi: ' num2str(mean(hexa_data_an.port_vis(qq).iti)) ' +/- ' num2str(std(hexa_data_an.port_vis(qq).iti))]);
    hexa_data_an.port_vis(qq).histX = 0.1:0.1:4;
    hexa_data_an.port_vis(qq).histY = hist(log10(hexa_data_an.port_vis(qq).iti),hexa_data_an.port_vis(qq).histX);
    plot(hexa_data_an.port_vis(qq).histX , cumsum(hexa_data_an.port_vis(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','2','3','4','5','6'},'Location','northwest');
        xlabel('Log10 Intervisit Interval'); ylabel('Count'); box off;
    end
    
    hexa_data_an.visits(qq,all_vis_inds(all_vis_ports==qq)) = 1;
    hexa_data_an.vi.avg(qq) = mean(hexa_data_an.port_vis(qq).iti);
    hexa_data_an.vi.std(qq) = std(hexa_data_an.port_vis(qq).iti);

end

hexa_data_an.all_rew_inds=all_rew_inds;
hexa_data_an.all_rew_ports=all_rew_ports;
hexa_data_an.all_vis_inds=all_vis_inds;
hexa_data_an.all_vis_ports=all_vis_ports;

subplot(133);
plot([1 max(log10(hexa_data_an.rw.avg))+1],[1 max(log10(hexa_data_an.rw.avg))+1],'k--'); hold on;
scatter(log10(hexa_data_an.vi.avg),log10(hexa_data_an.rw.avg),50,'filled');
axis([1 max(log10(hexa_data_an.rw.avg))+1 1 max(log10(hexa_data_an.rw.avg))+1]);
ylabel('Inter Reward Interval');
xlabel('Inter Visit Interval');
box off;
