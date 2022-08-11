%simulation of stem cell clonal expansion
tic
clear variables
clc


timepoint = [9000 15000 21000 29000 40000 52000 64000 103000];
x_matrix = importdata('timeseries_ep0.05_K420_N1000.mat');
load('variables.mat')
pickup = zeros(num_of_clones,length(timepoint));
for i = 1:length(timepoint)
    pickup(:,i) = x_matrix(:,timepoint(i))/n_openniche;
end
%%% pickup corresponds expdata


pickup_log = zeros(size(pickup));
for i = 1:length(pickup(:,1))
    for j = 1:length(pickup(1,:))
        if pickup(i,j) ~= 0
            pickup_log(i,j) = log(pickup(i,j));
        else
            pickup_log(i,j) = NaN;
        end
    end
end


%%% size heatmap
[pickup_max, max_ind] = max(pickup,[],2);
pickup_log_sort = sortrows(cat(2, max_ind, pickup_log, pickup_max));
pickup_log_sort_nz = zeros(nnz(pickup_log_sort(:,end)), length(pickup_log(1,:)));
j = 1;
for i = 1:length(pickup(:,1))
    if pickup_log_sort(i,end) ~= 0
        pickup_log_sort_nz(j,:) = pickup_log_sort(i,2:end-1);
        j = j+1;
    end
end
figure
imagesc(pickup_log_sort_nz,[-8,-2])
colorbar
saveas(gcf,'Fig5B.png')


%%% clone size distribution
figure
hold on
for i=1:8
    [n,xout]=hist(pickup(:,i),[0:0.005:0.06]);
    plot(xout,log(n),'- o')
    legend
end


%%% analytical solution of clone size distribution
r_plus = zeros(1,n_openniche+1);
r_minus = zeros(1,n_openniche+1);
for k = (0:n_openniche)
    r_plus(k+1) = (epsilon+lambda*k)*(n_openniche-k)/n_openniche/(num_of_clones*epsilon+lambda*n_openniche);
    r_minus(k+1) = (epsilon*(num_of_clones-1)+lambda*(n_openniche-k))*k/n_openniche/(num_of_clones*epsilon+lambda*n_openniche);
end
analytical_distribution = zeros(1,n_openniche+1);
h = zeros(1,n_openniche+1);
h(1) = 1;
for a = 2:n_openniche+1
    h(a) = h(a-1)*r_plus(a-1)/r_minus(a);
end
analytical_distribution(1) = 1/sum(h);
for a = 2:n_openniche+1
    analytical_distribution(a) = analytical_distribution(a-1)*r_plus(a-1)/r_minus(a);
end
x_axis = zeros(1,n_openniche+1);
for i = 1:n_openniche+1
    x_axis(i) = (i-1)/n_openniche;
end

plot(x_axis,log(num_of_clones*analytical_distribution*0.005*n_openniche),'Color','r')
%The coefficients are intended to compensate for the effect of bin width
%Δ=0.005*1000=5
%P_all(n)=K P(n)
%logP_all(Δ*i>n>=Δ*(i-1))= logΣ_n=Δ*(i-1)~ Δ*i Σ_k P_k(n)
%log P_all(n)=log {ΔK P(n)}

hold off
xlim([0 0.06])
ylim([-1 7])
saveas(gcf,'Fig5D.png')
toc