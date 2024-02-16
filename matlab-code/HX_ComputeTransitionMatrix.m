function [trans_mat] = HX_ComputeTransitionMatrix(tmp_visit_list,plotFlag,Nback)

visit_list = tmp_visit_list(1:end-Nback);

trans_mat = zeros(6,6);
for zz=1:6
    this_port_visits = find(visit_list==zz);
    % min(diff(this_port_visits))
    for mm=1:6
        trans_mat(zz,mm) = numel(find(tmp_visit_list(this_port_visits+Nback)==mm));
    end
end


% convert to probabilities
trans_mat = trans_mat ./ sum(sum(trans_mat));

if plotFlag>0
    figure(plotFlag); clf;
    exag = TNC_CreateRBColormap(8,'exag');
    imagesc(trans_mat,[0 0.2]); colormap(exag); axis equal; box off; colorbar;
end

