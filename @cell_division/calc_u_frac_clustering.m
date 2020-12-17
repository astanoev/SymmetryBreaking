function [fr_sst, clus_radius, mi, u_frac, range] = calc_u_frac_clustering(obj)
    fr_sst = zeros(length(obj.mss),3);
    clus_radius = zeros(length(obj.mss),1);
    mi = zeros(length(obj.mss),1);
    u_frac = nan(length(obj.mss),2000);
    for k=1:length(obj.mss)
        [~, ssl] = obj.mss(k).get_steady_state_labels();
        m = size(ssl,1); n = size(ssl,2);
        [X,Y] = meshgrid(1:n,1:m);
        coords = [X(:),Y(:)];
        states = ssl(:);
        
        dists_u_u = pdist(coords(states==1,:))';
        dists_u_v = pdist2(coords(states==1,:), coords(states>1,:));
        if isempty(dists_u_v); continue; end

        [X2,Y2] = meshgrid((0:n).^2,(0:m).^2);
        range = unique(sqrt(X2+Y2))';
        d = diff(range)/2; % set bin edges for histogram
        edges = [range(1)-d(1), range(1:end-1)+d, range(end)+d(end)];
        within_u = cumsum(histcounts(dists_u_u(:), edges));
        within_u = 2*within_u + nnz(states==1);
        within_v = cumsum(histcounts(dists_u_v(:), edges));
        u_frac(k,1:length(range)) = within_u./(within_u + within_v);
        
        fr_sst(k, :) = histcounts(states, 0.5:1:3.5)./length(states);
        fr_mid = (1+fr_sst(k,1))/2;
        ix = find(u_frac(k,:)<=fr_mid,1);
        clus_radius(k) = interp1(u_frac(k,[ix-1,ix]),range([ix-1,ix]),fr_mid);
        
        mi(k) = Morgan_I(obj.mss(k).neigh, states);
    end
    u_frac(:, (length(range)+1):end) = [];
end

function mi = Morgan_I(neigh, states)
    w = neighbourhood.hop_1_grid(neigh.m, neigh.n).A;
    w(eye(size(w))==1) = 0;
    fr_sst = nnz(states==1)/length(states);
    x = (states == 1)-fr_sst;
    num = x'*w*x;
    den = x'*x;
    mi = length(states)*num/(sum(sum(w))*den);
end

%u_nodes = states==1;
%v_nodes = states>1;
%n_u_states = nnz(u_nodes);
%dists = squareform(pdist(coords));
%dists_u_u = dists(u_nodes,u_nodes);
%dists_u_u = dists_u_u(:);
%dists_u_v = dists(u_nodes,v_nodes);
%range = unique(dists)';

%p1 = plot(ax,range, u_frac, 'linewidth', 2, 'DisplayName', num2str(k)); hold on;
%p2 = plot(ax,[20, range(end)], [fr_sst(k), fr_sst(k)], 'o', 'color', p1.Color);
%set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% legend();
% xlim([0,20]);
% ylim([0,1]);
% 
% figure;
% plot(1:length(obj.mss), clus, 'linewidth', 2);
% 
% figure;
% plot(1:length(obj.mss), fr_sst, 'linewidth', 2);

% u_frac = zeros(1, length(range));
% u_tot = zeros(1, length(range));
% u_frac_single = zeros(length(range), n_u_states);
% for j=1:length(range)
%     r = range(j);
%     within_u_single = sum(dists_u_u<=r, 2);
%     within_v_single = sum(dists_u_v<=r, 2);
%     weights = 1./(within_u_single +within_v_single);
%     weights = weights./sum(weights);
%     weights = 1;
%     within_u = sum(weights.*within_u_single);
%     within_v = sum(weights.*within_v_single);
%     %within_u = nnz(dists_u_u(:)<=r);%n_u_states+2*nnz(dists_u_u<=r);
%     %within_v = nnz(dists_u_v(:)<=r);
%     u_frac(j) = within_u/(within_u+within_v);
%     u_frac_single(j,:) = within_u_single./(within_u_single+within_v_single);
%     %u_frac(j) = sum(weights.*(within_u_single./(within_u_single+within_v_single)));
%     u_tot(j) = within_u+within_v;
% end
% %     figure; hold on;
% %     for i=1:n_u_states
% %         plot(range, u_frac_single(:,i)','DisplayName',strcat(num2str(coords(u_nodes(i),1)),',',num2str(coords(u_nodes(i),2))));
% %     end
% %     title(strcat(num2str(m),'x',num2str(n)));



% 
% figure;
% hold on;
% u_nodes = find(states==1);
% u_frac_range = zeros(length(u_nodes), length(range));
% parfor i=1:length(u_nodes)
%     u_frac_i = zeros(1, length(range));
%     for j=1:length(range)
%         r = range(j);
%         within = find(dists(u_nodes(i),:)<=r);
%         %u_frac_range(i,j) = nnz(states(within)==1)/length(within);
%         u_frac_i(j) = nnz(states(within)==1)/length(within);
%     end
%     u_frac_range(i,:) = u_frac_i;
%     %plot(range, u_frac_range(i,:));
% end
% %plot(range, u_frac_range);
% plot(range, mean(u_frac_range, 1), 'linewidth', 2);
% plot(range, median(u_frac_range, 1), 'linewidth', 2);
% u_frac_median = median(u_frac_range, 1);
% fr_sst = nnz(states==1)/length(states)
% mid_fr = (1+fr_sst)/2;
% ix = find(u_frac_median<=mid_fr,1);
% range(ix)