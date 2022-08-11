tic
clear variables
clc

load('variables.mat')
x_matrix = importdata('timeseries_ep0.9_K10_N100.mat');
epsilon = 0.9;
lambda = 1;


tmax = 500;
figure
for k=1:num_of_clones
    subplot(num_of_clones,1,k)
    plot(x_matrix(k,:))
    xlim([0 tmax*(lambda*n_openniche+epsilon*num_of_clones)])
    ylim([0,n_openniche])
    title(['clone ',num2str(k)])
    set(gca,'xtick',[],'ytick',[])
end


toc