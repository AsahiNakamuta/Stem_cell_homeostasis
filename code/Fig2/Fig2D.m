tic
clear variables
clc

load("variables.mat")

%lists of proliferation rate
epsilon_list = [0 0.01 0.03 0.05 1];
lambda_list  = [1 1 1 1 0];

%number of simulation trials for each parameter set
m = 1000;
%number of clones (K in manuscript)
num_of_clones=10;
%number of competitive stem cells in the open layer (N in manuscript)
n_openniche = 100;
%list of number of iteration
iter_list = 500*(n_openniche*lambda_list + num_of_clones*epsilon_list);

%time-series data of the clone #1
x_clone1 = importdata("timeseries_monoclonalconversion.mat");

%time-series data of the probability of monoclonal conversion
mono = zeros(length(epsilon_list),max(iter_list));

for i_para = 1:length(epsilon_list)
    for i = 1:iter_list(i_para)
        for j = 1:m
            if x_clone1(i_para,j,i)==n_openniche
                mono(i_para,i) = mono(i_para,i)+1;
            end
        end
        mono(i_para,i) = mono(i_para,i)/nnz(x_clone1(i_para,:,i));
    end
end

hold on
for e = 1:length(epsilon_list)
    plot(mono(e,1:iter_list(e)))
end
legend('ε=0','ε=0.01','ε=0.1','ε=1','λ=0')
hold off
toc