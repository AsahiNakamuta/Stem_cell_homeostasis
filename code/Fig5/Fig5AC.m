clear variables
clc
tic


expdata = readmatrix('mmc6.xlsx','Sheet','2RC003','Range','Z2:AH615');


expdata_log = zeros(size(expdata));
for i = 1:length(expdata(:,1))
    for j = 1:length(expdata(1,:))
        if expdata(i,j) ~= 0
            expdata_log(i,j) = log(expdata(i,j));
        else
            expdata_log(i,j) = NaN;
        end
    end
end


%%% size heatmap
[expdata_max, max_ind] = max(expdata,[],2);
expdata_log_sort = sortrows(cat(2, max_ind, expdata_log));
expdata_max = sortrows(cat(2, max_ind, expdata_max));
expdata_log_sort_nz = zeros(nnz(expdata_max(:,2)), length(expdata_log_sort(1,:)));
j = 1;
for i = 1:length(expdata_max(:,2))
    if expdata_max(i,2) ~= 0
        expdata_log_sort_nz(j,:) = expdata_log_sort(i,:);
        j = j+1;
    end
end
figure
imagesc(expdata_log_sort_nz(:,3:end),[-8,-2])
colorbar
saveas(gcf,'Fig5A.png')


%%% clone size distribution
figure
hold on
for i=2:9
    [n,xout]=hist(expdata(:,i),[0:0.005:0.06]);
    plot(xout,log(n),'- o')
    legend
end
hold off
ylim([-1 7])
saveas(gcf,'Fig5C.png')

toc