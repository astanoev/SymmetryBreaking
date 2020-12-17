neigh_cl = @neighbourhood.range_grid;
neigh_range = 1;
%neigh_args = {};
model_cl = @models.u_v_dn_model;
alfa1s = 2.00:0.005:2.4;%1.8:.005:2.8;
Ns = 2.^[1:8];
%Ns = 1:128;
n_reps = 200;

par_comb = utils.combvec(alfa1s, Ns, 1:n_reps);
out_par = zeros(size(par_comb,2),1);
out = zeros(length(alfa1s),length(Ns),n_reps);
first_unst = inf(length(alfa1s),n_reps);
first_unst_2 = inf(length(alfa1s),n_reps);

for i=1:length(alfa1s)
    alfa1 = alfa1s(i);
    disp(alfa1);
    par = struct('alfa1',alfa1);
    model = model_cl(neighbourhood.self_grid(1,1),'par',par);
    %model.update_parameters(par);
    mlp = model.mlp;
    
    for rep=1:n_reps
        flag = 0;
        Ns_rep = Ns;
        while true
            for j=1:length(Ns_rep)
                N = Ns_rep(j);
                [m,n] = neighbourhood.neighbourhood.get_m_n_N(N);
            
                neigh = neigh_cl(m,n,neigh_range);
                model.neigh = neigh;
                %model = model_cl(neigh,'mlp',mlp,'par',par,'mock',true);
                %model.update_parameters(par);
                %model.set_initial_cond(mlp);
                [~, stable] = model.jacobian(model.mlp_mat(:));
                if stable==0; flag = flag+1; break; end
            end
            if flag==1 && isinf(first_unst_2(i,rep))
                first_unst_2(i,rep) = N;
                first_unst(i,rep) = N;
                if j==1
                    flag = flag+1;
                else
                    Ns_rep = (Ns(j-1)+1):(N-1); % for first_unst search
                end
                %Ns = 1:(N-1);
                break;
            elseif flag==2
                first_unst(i,rep) = N;
                %break;
            end
            if flag~=1 || first_unst(i,rep)~=N || isempty(Ns_rep); break; end
        end
    end
end

% for p_i=1:size(par_comb,2)
%     params = num2cell(par_comb(:,p_i));
%     [i,j,rep] = params{:};
%     out(i,j,rep) = out_par(p_i);
% end

% out_mean = mean(out,3);
% figure;
% imshow(out_mean');

figure; hold on; box on;
alfa1s_maxprob = zeros(length(alfa1s),2);
for i=1:length(alfa1s)
    alfa1 = alfa1s(i);
    fu2 = first_unst_2(i,:);
    [fu2_un, ~, ic] = unique(fu2);
    fu2_un(isinf(fu2_un)) = 70;
    a_counts = accumarray(ic,1);
    ix = find(cumsum(a_counts)>=0.5*n_reps,1,'first');
    %[~, ix] = max(a_counts);
    alfa1s_maxprob(i,:) = [alfa1, fu2_un(ix)];
    scatter(alfa1*ones(size(fu2_un)),fu2_un,0.5*a_counts,'b','filled');
end
plot(alfa1s_maxprob(:,1),alfa1s_maxprob(:,2),'b-');
% plot(alfa1s,mean(first_unst_2,2),'b-');
% plot(alfa1s,first_unst,'ro');
% plot(alfa1s,mean(first_unst,2),'r-');
xlim([alfa1s(1),alfa1s(end)]);
ylim([0,max(first_unst_2(~isinf(first_unst_2)))]);