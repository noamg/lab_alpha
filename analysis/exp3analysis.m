folder = 'C:\Users\Dan Wexler\Documents\school\year 3\lab3\alpha decay\alpha.git\exp3\';
files = dir([folder,'*.Spe']);

for ind = 1:length(files)
    unsortInd(ind) = str2num(files(ind).name(1:end-4));
end
[~,sortInd] = sort(unsortInd);

for ind = 1:length(files)
    exp3.time(ind,1:2) = dlmread([folder,files(sortInd(ind)).name],' ',[9,0,9,1]);
    exp3.data(ind,1:2048) = dlmread([folder,files(sortInd(ind)).name],'\t',[12,0,12+2047,0]);
    exp3.data(ind,1:100) = 0;
    %data: matrix: lines are according to times, cols according to channal,
    %values are counts
end

%diff the counts by time, divide by dt to get rate
rate = bsxfun(@rdivide, diff(exp3.data, 1, 1), diff(exp3.time(:,1)));
exp3error = bsxfun(@rdivide, sqrt(diff(exp3.data, 1, 1)), diff(exp3.time(:,1)));
t_means = 0.5 * (exp3.time(1:end-1,1) + exp3.time(2:end, 1));
%find peaks, peaks area
peak_l = [906, 940, 990, 1313];
areas = zeros(length(t_means), 4);

for line = 1:length(t_means)

    temp_fit = skewGaussFit_4(rate(line, :), peak_l);
    areas(line, 1) = temp_fit.a1;
    areas(line, 2) = temp_fit.a2;
    areas(line, 3) = temp_fit.a3;
    areas(line, 4) = temp_fit.a4;
%     tempErr = diff(confint(temp_fit));
%     areaErr(line,:) = tempErr(1:4)/2;
end

areasErrRelative = exp3error(:,peak_l)./rate(:,peak_l);

figure
errorbar(repmat(t_means,1,4), areas,areasErrRelative.*areas)
%plot peaks by time (4 series)
figure
plot(t_means, areas)
% errorbar(repmat(t_means,1,4), areas,areaErr)
figure
v = areas(:,3);
plot(t_means, log(v))

foo = @(a,tau, x) a*exp(-x./tau);

fit(t_means, v, foo)


ix = [1,10,20];
figure
for ind = 1:length(ix)
    subplot(3,1,ind)
    hold on
    h = plot(bin2E(0:2047),rate(ix(ind),:));
    axis([0,bin2E(2050),0,2])
    grid
    set(h,'displayname',['t = ',num2str(round(t_means(ix(ind))/60,0),2),' min'])
    legend('show')
%     plot(bin2E(0:2047),rate(ix(ind),:)-exp3error(ix(ind),:),'--')
%     plot(bin2E(0:2047),rate(ix(ind),:)+exp3error(ix(ind),:),'--')
end





