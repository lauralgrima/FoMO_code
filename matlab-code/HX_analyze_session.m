function [hexa_data_an] = HX_analyze_session(hexa_data,session)

all_rew_inds = find(hexa_data.rewarded==1 & hexa_data.session_n==session);
all_rew_ports = hexa_data.port_n(all_rew_inds);

all_vis_inds = find(hexa_data.unique_vis==1 & hexa_data.session_n==session);
all_vis_ports = hexa_data.port_n(all_vis_inds);

figure(57); clf;
cat_map = TNC_CreateRBColormap(8,'mapb');

for qq=unique(all_rew_ports)

    subplot(121)
    hexa_data_an.port_rew(qq).ts = all_rew_inds(all_rew_ports==qq);
    hexa_data_an.port_rew(qq).iti = diff(hexa_data_an.port_rew(qq).ts);
    disp(['Port ' num2str(qq) ' mean iti: ' num2str(mean(hexa_data_an.port_rew(qq).iti)) ' +/- ' num2str(std(hexa_data_an.port_rew(qq).iti))]);
    hexa_data_an.port_rew(qq).histX = 0.1:0.1:4;
    hexa_data_an.port_rew(qq).histY = hist(log10(hexa_data_an.port_rew(qq).iti),hexa_data_an.port_rew(qq).histX);
    plot(hexa_data_an.port_rew(qq).histX , cumsum(hexa_data_an.port_rew(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','2','3','4','5','6'},'Location','northwest');
        xlabel('Log10 Interreward Interval'); ylabel('Count'); box off;
    end

    subplot(122)
    hexa_data_an.port_vis(qq).ts = all_vis_inds(all_vis_ports==qq);
    hexa_data_an.port_vis(qq).iti = diff(hexa_data_an.port_vis(qq).ts);
    disp(['Port ' num2str(qq) ' mean iti: ' num2str(mean(hexa_data_an.port_vis(qq).iti)) ' +/- ' num2str(std(hexa_data_an.port_vis(qq).iti))]);
    hexa_data_an.port_vis(qq).histX = 0.1:0.1:4;
    hexa_data_an.port_vis(qq).histY = hist(log10(hexa_data_an.port_vis(qq).iti),hexa_data_an.port_vis(qq).histX);
    plot(hexa_data_an.port_vis(qq).histX , cumsum(hexa_data_an.port_vis(qq).histY) , 'color' , cat_map(qq,:) , 'linewidth', 2); hold on;
    if qq==max(all_rew_ports)
        legend({'1','2','3','4','5','6'},'Location','northwest');
        xlabel('Log10 Intervisit Interval'); ylabel('Count'); box off;
    end
    
end


hexa_data_an.port_chk
hexa_data_an.visits