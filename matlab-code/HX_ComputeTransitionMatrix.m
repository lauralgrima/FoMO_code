function [trans_mat] = HX_ComputeTransitionMatrix(visit_list,plotFlag)

trans_mat = zeros(6,6);
for zz=1:6
    this_port_visits = find(visit_list(2:end)==zz);
    trans_mat(zz,:) = histcounts(visit_list(this_port_visits),6);
end

if plotFlag>0
    figure(plotFlag); clf;
    exag = TNC_CreateRBColormap(8,'exag');
    imagesc(trans_mat./numel(visit_list),[0 0.2]); colormap(exag); axis equal; box off;
end