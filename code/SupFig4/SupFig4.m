%scaling law in clone size

tic
clear variables
clc

load('variables.mat')
x_clone1 = importdata("clone1timeseries_ep10_K10_N100.mat");

epsilon = 10;
lambda = 1;


temp = zeros(3,n_openniche);
unscale = zeros(3,n_openniche);
scale = zeros(3,n_openniche);
x_axis = zeros(3,n_openniche);
for t = 1:3
    pickup = round(2*t*(epsilon*num_of_clones+lambda*n_openniche));
    for i = 1:m
        if x_clone1(i,pickup)~=0
            temp(t,x_clone1(i,pickup)) = temp(t,x_clone1(i,pickup))+1;
        end
    end
    for j=1:n_openniche
        unscale(t,j) = temp(t,j)/nnz(x_clone1(:,pickup));
        scale(t,j) = temp(t,j)/nnz(x_clone1(:,pickup))*sum(x_clone1(:,pickup))/nnz(x_clone1(:,pickup));
        x_axis(t,j) = j/(sum(x_clone1(:,pickup))/nnz(x_clone1(:,pickup)));
    end
end

figure
semilogy(0,0)
hold on
for t = 1:3
    semilogy(unscale(t,:))
end
hold off
%xlim([0,100])
%ylim([0.00001, 0.1])
figure
semilogy(0,0)
hold on
for t = 1:3
    semilogy(x_axis(t,:),scale(t,:))
end
hold off
%xlim([0,9])

toc