clear all
close all

M = readmatrix('cd_data.xlsx');

[B,I] = sort(M);
I = I';
trial_num = B(:,1)';
cue_DA = M(I(1,:),2)';
rew_DA = M(I(1,:),3)';
total_DA=cue_DA+rew_DA;

avg_cue_DA = movmean(cue_DA,10);
avg_rew_DA = movmean(rew_DA,10);
avg_rewcue_DA = movmean(total_DA,10);

figure;
subplot(1,3,1)
scatter(trial_num,cue_DA)
hold on
plot(trial_num,avg_cue_DA)
hold on
plot(trial_num,0*trial_num, 'k--')
hold off
xlabel('Training trial')
ylabel('change in DA at cue (Hz)')
subplot(1,3,2)
scatter(trial_num,rew_DA)
hold on
plot(trial_num,avg_rew_DA)
hold on
plot(trial_num,0*trial_num, 'k--')
hold off
xlabel('Training trial')
ylabel('change in DA at reward (Hz)')
subplot(1,3,3)
scatter(trial_num,total_DA)
hold on
plot(trial_num,avg_rewcue_DA)
plot(trial_num,0*trial_num, 'k--')

hold off
xlabel('Training trial')
ylabel('change in DA at reward and cue combined (Hz)')

figure(2)
scatter(trial_num,total_DA,40,'filled')
hold on
plot(trial_num,avg_rewcue_DA,'k','LineWidth',2)
plot(trial_num,0*trial_num, 'k--')


xlabel('Training trial')
ylabel('Total DA change (cue+reward) (Hz)')

basal_DA=total_DA(1:10);
P=zeros([1,9]);
H=zeros([1,9]);
for i=1:8
    test_DA=total_DA(1+i*10:(i+1)*10);
    [P(i),H(i)]=ranksum(basal_DA,test_DA);
    if H(i)==1
        xx=mean(trial_num(1+i*10:(i+1)*10));
        plot(xx,23,'*','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',1);
    end %if
end
test_DA=total_DA(87:96);
[P(9),H(9)]=ranksum(basal_DA,test_DA);

if H(9)==1
    xx=mean(trial_num(87:96));
    plot(xx,23,'*','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',1);
end

    
% for i=1:11
%     test_DA=total_DA(1+i*8:(i+1)*8);
%     [P(i),H(i)]=ranksum(basal_DA,test_DA);
%     if H(i)==1
%         xx=mean(trial_num(1+i*8:(i+1)*8));
%         plot(xx,22,'*','MarkerEdgeColor','k','MarkerSize',10,'LineWidth',1);
%     end %if
% end