tic
clear variables
clc


load('variables.mat')
num_of_clones = 10;
n_openniche = 100;
lambda = 1;

sizedist = zeros(3,n_openniche+1);
burstprob = zeros(3,n_openniche);
duration = zeros(3,n_openniche-1);


for i = 1:3
    epsilon = i*0.4 - 0.3;
    r_plus = zeros(1,n_openniche+1);
    r_minus = zeros(1,n_openniche+1);
    for k = 0:n_openniche
        r_plus(k+1) = (epsilon+lambda*k)*(n_openniche-k)/n_openniche/(num_of_clones*epsilon+lambda*n_openniche);
        r_minus(k+1) = (epsilon*(num_of_clones-1)+lambda*(n_openniche-k))*k/n_openniche/(num_of_clones*epsilon+lambda*n_openniche);
    end


    %clone size distribution
    hoge = zeros(1,n_openniche+1);
    hoge(1) = 1;
    for a = 2:n_openniche+1
        hoge(a) = hoge(a-1)*r_plus(a-1)/r_minus(a);
    end
    sizedist(i,1) = 1/sum(hoge);
    for a = 2:n_openniche+1
        sizedist(i,a) = sizedist(i,a-1)*r_plus(a-1)/r_minus(a);
    end


    %Probability of burst generation of each height
    %to calculate 'burstprob', calculate 'temp1' and 'temp2' first.
    temp1 = zeros(1,n_openniche);
    temp2 = zeros(n_openniche);
    
    temp1(1) = 1;
    for j = 2:n_openniche
        temp1(j) = temp1(j-1)*r_minus(j)/r_plus(j);
    end
    for j = 1:n_openniche
        temp2(j,j) = 1;
        for k = j+1:n_openniche
            temp2(j,k) = temp2(j,k-1)*r_plus(k)/r_minus(k);
        end
    end

    %then 'burstprob'
    for j = 1:n_openniche
        burstprob(i,j) = r_plus(1)/sum(temp1(1:j))/sum(temp2(:,j));
    end


    %Duration of burst generation of each height
    %Duration of the first half of a burst
    duration_forward = zeros(1,n_openniche);
    for j = 2:n_openniche
        syms var_temp
        temp3 = sym('temp3',[1 j]);
        temp3(1) = var_temp;
        for k = 2:j
            temp3(k) = (r_minus(k)*temp3(k-1) - sum(temp1(1:k-1))/sum(temp1(1:j)))/r_plus(k);
        end
        eqn = sum(temp3)==0;
        answer = solve(eqn,var_temp);
        duration_forward(j) = double(answer*sum(temp1(1:j))+1);
    end
    %Duration of the last half of a burst
    duration_backward = zeros(1,n_openniche);
    for j = 2:n_openniche
        syms var_temp
        temp4 = sym('temp4',[1 j]);
        temp4(j) = -var_temp;
        for k = j-1:-1:1
            temp4(k) = (r_plus(k+1)*temp4(k+1) + sum(temp2(k+1:j,j))/sum(temp2(:,j)))/r_minus(k+1);
        end
        eqn = sum(temp4)==0;
        answer = solve(eqn,var_temp);
        duration_backward(j-1) = double(answer*sum(temp2(:,j))+1);
    end
    %total duration
    duration(i,:) = duration_forward(1:n_openniche-1) + duration_backward(1:n_openniche-1);
end


figure
hold on 
for i = 1:3
    plot(log(burstprob(i,:)))
end
hold off
figure
hold on 
for i = 1:3
    plot(duration(i,:))
end
hold off


%save
filename = ['_K', num2str(num_of_clones), '_N', num2str(n_openniche), '.mat'];
save(['ana_sizedist',filename],'sizedist')
save(['ana_burstprob',filename],'burstprob')
save(['ana_duration',filename],'duration')
toc