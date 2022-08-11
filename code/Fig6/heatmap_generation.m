%%% heatmap of the Shannon index

tic
clear variables
clc

%initial conditions and parameters

%number of iteration
iter=100000;
%number of clones (K in manuscript)
num_of_clones=420;
%matrix containing time series data
x_matrix=zeros(num_of_clones,iter);
%proliferation rate of competitive stem cells
lambda=1;
%number of rows or columns of heatmap
n_heatmap = 100;
% matrix
shannon_heatmap = zeros(n_heatmap,n_heatmap);
duration_heatmap = zeros(n_heatmap,n_heatmap);

%number of competitive stem cells in the open layer (N in manuscript)
n_openniche = 100;
%initional condition
init_size = floor(n_openniche/num_of_clones);
x_matrix(:,1) = init_size*ones(num_of_clones,1);
x_matrix(1:n_openniche-init_size*num_of_clones,1) = x_matrix(1:n_openniche-init_size*num_of_clones,1) + 1;
%proliferation rate of matster stem cells
epsilon=0.3;
for i_openniche = 1:n_heatmap
    for i_epsilon = 1:n_heatmap
        %proliferation rate of matster stem cells
        epsilon=i_epsilon/n_heatmap;
        %number of cells in the open niche
        n_openniche = i_openniche/n_heatmap*1000;
        %initional condition
        init_size = floor(n_openniche/num_of_clones);
        x_matrix(:,1) = init_size*ones(num_of_clones,1);
        x_matrix(1:n_openniche-init_size*num_of_clones,1) = x_matrix(1:n_openniche-init_size*num_of_clones,1) + 1;
        
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

        %%%shannon index
        shannon = 0;
        for i = 1:num_of_clones
            for j = 1:iter
                if x_matrix(i,j)~=0
                   shannon = shannon - x_matrix(i,j)/n_openniche*log(x_matrix(i,j)/n_openniche);
                end
            end
        end
        shannon_heatmap(i_openniche,i_epsilon) = shannon/iter;

        %%%average duration
        nburst = 0;
        duration = 0;
        for i = 1:num_of_clones
            j = 1;
            while j < iter
                if x_matrix(i,j) == 0
                    j = j+1;
                else
                    start = j;
                    j = j+1;
                    while j < iter
                        if x_matrix(i,j) ~= 0
                            j = j+1;
                        else
                            nburst = nburst+1;
                            stop = j;
                            duration = duration + stop - start;
                            j = j+1;
                            break
                        end
                    end
                end
            end
        end
        duration_heatmap(i_openniche,i_epsilon) = duration/nburst/(lambda*n_openniche+epsilon*num_of_clones);

    end
end
save('shannon_heatmap.mat', "shannon_heatmap")
save('duration_heatmap.mat', "duration_heatmap")
x = [10 1000];
y = [0.01 1];
figure
imagesc(x,y,shannon_heatmap)
colorbar
figure
imagesc(x,y,duration_heatmap)
colorbar

toc