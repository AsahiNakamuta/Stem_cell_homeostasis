%Fig3A

tic
clear variables
clc

simu_sizedist = importdata('simu_sizedist_K10_N100.mat');
simu_forward_burstprob = importdata('simu_forward_burstprob_K10_N100.mat');
simu_backward_burstprob = importdata('simu_backward_burstprob_K10_N100.mat');
simu_forward_duration = importdata('simu_forward_duration_K10_N100.mat');
simu_backward_duration = importdata('simu_backward_duration_K10_N100.mat');
ana_sizedist = importdata('ana_sizedist_K10_N100.mat');
ana_burstprob = importdata("ana_burstprob_K10_N100.mat");
ana_duration = importdata('ana_duration_K10_N100.mat');
load('variables.mat')
x_axis = [9 19 29 39 49 59 69 79 89 99];

figure
hold on
for i = 1:3
    plot(log(ana_sizedist(i,:)))
    plot(log(simu_sizedist(i,:)),'o','MarkerSize',3)
end
xlim([0,100])
ylim([-30,0])
hold off

figure
hold on
for i = 1:3
    plot(log(ana_burstprob(i,:)))
    plot(x_axis,log(simu_forward_burstprob(i,:).*simu_backward_burstprob(i,:)),'o','MarkerSize',5)
end
xlim([0,100])
ylim([-30,0])
hold off

figure
hold on
for i = 1:3
    plot(ana_duration(i,:)/(num_of_clones*(i*0.4-0.3) + n_openniche*lambda));
    plot(x_axis,(simu_forward_duration(i,:)+simu_backward_duration(i,:))/(num_of_clones*(i*0.4-0.3) + n_openniche*lambda),'o','MarkerSize',5)
end
xlim([0,100])
%ylim([0,250])
hold off