%simulation of stem cell clonal expansion
tic
clear variables
clc


%initial conditions and parameters

%number of iteration
iter=103000;
%number of clones (K in manuscript)
num_of_clones=10;
%number of competitive stem cells in the open layer (N in manuscript)
n_openniche = 100;
%matrix containing time series data
x_matrix=zeros(num_of_clones,iter);
%initional condition
init_size = floor(n_openniche/num_of_clones);
x_matrix(:,1) = init_size*ones(num_of_clones,1);
x_matrix(1:n_openniche-init_size*num_of_clones,1) = x_matrix(1:n_openniche-init_size*num_of_clones,1) + 1;
%proliferation rate of matster stem cells
epsilon=0;
%proliferation rate of competitive stem cells
lambda=1;


%time evolution
for i=1:iter-1
    
    x_temp=x_matrix(:,i);
%lists of the probabilities that each clone will be selected for proliferation/differentiation.
    f_pro_list=(lambda*x_temp+epsilon)/(lambda*n_openniche+num_of_clones*epsilon);
    f_diff_list=x_temp/n_openniche;
%choose one clone to proliferate
    threshold_pro_list=zeros(num_of_clones-1,1);
    threshold_pro_list(1)=f_pro_list(1);
    
    for k=2:num_of_clones-1
        threshold_pro_list(k)=threshold_pro_list(k-1)+f_pro_list(k);
    end
    
    rand_num1=rand;
    pro_type=sum(threshold_pro_list<=rand_num1)+1;
    
    x_temp(pro_type)=x_temp(pro_type)+1;
%choose one clone to differenciate
    threshold_diff_list=zeros(num_of_clones-1,1);
    threshold_diff_list(1)=f_diff_list(1);
    
    for k=2:num_of_clones-1
        threshold_diff_list(k)=threshold_diff_list(k-1)+f_diff_list(k);
    end
    
    rand_num2=rand;
    diff_type=sum(threshold_diff_list<=rand_num2)+1;
    
    x_temp(diff_type)=x_temp(diff_type)-1;
    
    x_matrix(:,i+1)=x_temp;
    
end


%%save
filename = ['timeseries', '_ep', num2str(epsilon), '_K', num2str(num_of_clones), '_N', num2str(n_openniche), '.mat'];
save(filename,'x_matrix')
save('variables.mat')


toc