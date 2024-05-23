clear all
close all

load('ff_values.mat');

value_list = [.25 .5 .75 1 1.25 1.5];

[~,~,stats] = anova2(ff_values',1,"off");
c1 = multcompare(stats);

tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


% d = randi(10, 20, 3);
% figure(1)
% scatter(value_list,ff_values)
% yt = get(gca, 'YTick');
% axis([xlim    0  ceil(max(yt)*1.2)])
% xt = get(gca, 'XTick');
% hold on
% plot(xt([2 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([2 3])), max(yt)*1.15, '*k')
% hold off