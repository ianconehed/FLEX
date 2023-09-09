%% average weights plotting
figure('Position',[800 100 600 600])
subplot(3,1,1)
plot(ff_vect(1,:),'x');
if num_columns == 2
    hold on
    plot(ff_vect1(1,:),'x');
end
xlabel('Trial Number');
ylabel('Mean FF Weight');
subplot(3,1,2)
plot(rec_vect(1,:),'x');
if num_columns == 2
    hold on
    plot(rec_vect1(1,:),'x');
end
xlabel('Trial Number');
ylabel('Mean rec Weight');
subplot(3,1,3)
plot(m_vect(1,:),'x');
if num_columns == 2
    hold on
    plot(m_vect1(1,:),'x');
end
xlabel('Trial Number');
ylabel('Mean inh Weight');
%% means for plotting

sample_len = 30;
trial_sel_1 = ceil(1); 
trial_sel_2 = ceil(num_trials/4 + 15);
trial_sel_3 = ceil(num_trials/2.6); 
trial_sel_4 = ceil(num_trials/2); 
trial_sel_5 = ceil(num_trials-lengthie-1); 


mean_VTA_plot1 = mean(R_it_l(N-num_VTA+1:N,:,trial_sel_1),1);
mean_VTA_plot2 = mean(R_it_l(N-num_VTA+1:N,:,trial_sel_2),1);
mean_VTA_plot3 = mean(R_it_l(N-num_VTA+1:N,:,trial_sel_3),1);
mean_VTA_plot4 = mean(R_it_l(N-num_VTA+1:N,:,trial_sel_4),1);
mean_VTA_plot5 = mean(R_it_l(N-num_VTA+1:N,:,trial_sel_5),1);

mean_R_it_plot1 = mean(R_it_l(:,:,trial_sel_1:trial_sel_1+sample_len),3);
mean_R_it_plot2 = mean(R_it_l(:,:,trial_sel_2-sample_len/2:trial_sel_2+sample_len/2),3);
mean_R_it_plot3 = mean(R_it_l(:,:,trial_sel_3:trial_sel_3+sample_len),3);
mean_R_it_plot4 = mean(R_it_l(:,:,trial_sel_4:trial_sel_4+sample_len),3);
mean_R_it_plot5 = mean(R_it_l(:,:,trial_sel_5:trial_sel_5+sample_len),3);




%% plot dopamine

% figure
% plot(dopa)
% hold on
% plot(dopa_ff)
% figure
% plot(sum(temp_t,1)/num_VTA)
% figure
% plot(dopa.*T_pt)
% hold on
% plot(dopa_ff.*T_dt)

%% average fr

figure('Position',[100 100 700 700])
subplot(5,1,1)
plot(1:dt:t_total,mean_VTA_plot1)
ylabel('Firing rate')
xlabel('Time(ms)')
title("Mean VTA neuron firing rate for trial "+ num2str(trial_sel_1))
xlim([0.01 t_total]);
ylim([0.01 max([mean_VTA_plot1 mean_VTA_plot2 mean_VTA_plot3 mean_VTA_plot4 mean_VTA_plot5]) + 5]);
subplot(5,1,2)
plot(1:dt:t_total,mean_VTA_plot2)
ylabel('Firing rate')
xlabel('Time(ms)')
title("Mean VTA neuron firing rate for trial "+ num2str(trial_sel_2))
xlim([0.01 t_total]);
ylim([0.01 max([mean_VTA_plot1 mean_VTA_plot2 mean_VTA_plot3 mean_VTA_plot4 mean_VTA_plot5]) + 5]);
subplot(5,1,3)
plot(1:dt:t_total,mean_VTA_plot3)
ylabel('Firing rate')
xlabel('Time(ms)')
title("Mean VTA neuron firing rate for trial "+ num2str(trial_sel_3))
xlim([0.01 t_total]);
ylim([0.01 max([mean_VTA_plot1 mean_VTA_plot2 mean_VTA_plot3 mean_VTA_plot4 mean_VTA_plot5]) + 5]);
subplot(5,1,4)
plot(1:dt:t_total,mean_VTA_plot4)
ylabel('Firing rate')
xlabel('Time(ms)')
title("Mean VTA neuron firing rate for trial "+ num2str(trial_sel_4))
xlim([0.01 t_total]);
ylim([0.01 max([mean_VTA_plot1 mean_VTA_plot2 mean_VTA_plot3 mean_VTA_plot4 mean_VTA_plot5]) + 5]);
subplot(5,1,5)
plot(1:dt:t_total,mean_VTA_plot5)
ylabel('Firing rate')
xlabel('Time(ms)')
title("Mean VTA neuron firing rate for trial "+ num2str(trial_sel_5))
xlim([0.01 t_total]);
ylim([0.01 max([mean_VTA_plot1 mean_VTA_plot2 mean_VTA_plot3 mean_VTA_plot4 mean_VTA_plot5]) + 5]);

% figure('Position',[300 300 1000 1000])
% subplot(2,1,1)
% LineFormat = struct();
% LineFormat.Color = 'yellow';
% plotSpikeRaster(sc_R_it(:,:,trial_sel_1) == 1,'LineFormat',LineFormat);
% set(gca,'color',[0 0 0])
% ylabel('Neuron number')
% xlabel('Time(ms)')
% title('spike raster  (200-300 are VTA neurons)')
% subplot(2,1,2)
% plot(1:dt:t_total,mean_VTA_plot2)
% ylabel('Firing rate')
% xlabel('Time(ms)')
% title('Mean VTA neuron firing rate, last trial')
% xlim([0.01 t_total]);
% ylim([0.01 max(mean_VTA_plot2) + 5]);


%% individual averaged fr

figure('Position',[50 300 1800 1000])
subplot(1,5,1)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_plot1(i,:))
end
view(0,85)
xlim([0 t_total])
% xlabel('Time(ms)')
% ylabel('Neuron number')
% zlabel('Firing rate (Hz)')
% title("VTA neurons, averaged over trials " + num2str(trial_sel_1) +...
%     " to " + num2str(trial_sel_1+sample_len))
colororder(white(100)-1)
hold off

subplot(1,5,2)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_plot2(i,:))
end
view(0,85)
xlim([0 t_total])
colororder(white(100)-1)
hold off

subplot(1,5,3)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_plot3(i,:))
end
view(0,85)
xlim([0 t_total])
colororder(white(100)-1)
hold off

subplot(1,5,4)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_plot4(i,:))
end
view(0,85)
xlim([0 t_total])
colororder(white(100)-1)
hold off

subplot(1,5,5)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_plot5(i,:))
end
view(0,85)
xlim([0 t_total])
colororder(white(100)-1)
hold off


%% mean plotting

max_vect = max([mean(mean_R_it_plot1(N-num_VTA+1:N,:)) mean(mean_R_it_plot2(N-num_VTA+1:N,:))...
     mean(mean_R_it_plot3(N-num_VTA+1:N,:)) mean(mean_R_it_plot4(N-num_VTA+1:N,:))...
      mean(mean_R_it_plot5(N-num_VTA+1:N,:))]);
figure
subplot(5,1,1)
plot(1:dt:t_total,mean(mean_R_it_plot1(N-num_VTA+1:N,:)))
title("Average over all VTA cells and trials " + num2str(trial_sel_1) +...
    " to " + num2str(trial_sel_1+sample_len))
xlabel('time')
ylabel('Firing rate')
ylim([0.01 max_vect + 5]);
subplot(5,1,2)
plot(1:dt:t_total,mean(mean_R_it_plot2(N-num_VTA+1:N,:)))
title("Average over all VTA cells and trials " + num2str(trial_sel_2) +...
    " to " + num2str(trial_sel_2+sample_len))
xlabel('time')
ylabel('Firing rate')
ylim([0.01 max_vect + 5]);
subplot(5,1,3)
plot(1:dt:t_total,mean(mean_R_it_plot3(N-num_VTA+1:N,:)))
title("Average over all VTA cells and trials " + num2str(trial_sel_3) +...
    " to " + num2str(trial_sel_3+sample_len))
xlabel('time')
ylabel('Firing rate')
ylim([0.01 max_vect + 5]);
subplot(5,1,4)
plot(1:dt:t_total,mean(mean_R_it_plot4(N-num_VTA+1:N,:)))
title("Average over all VTA cells and trials " + num2str(trial_sel_4) +...
    " to " + num2str(trial_sel_4+sample_len))
xlabel('time')
ylabel('Firing rate')
ylim([0.01 max_vect + 5]);
subplot(5,1,5)
plot(1:dt:t_total,mean(mean_R_it_plot5(N-num_VTA+1:N,:)))
title("Average over all VTA cells and trials " + num2str(trial_sel_5) +...
    " to " + num2str(trial_sel_5+sample_len))
xlabel('time')
ylabel('Firing rate')
ylim([0.01 max_vect + 5]);




%% AUC plotting
% 
% 
% baseline = hist_vect1(:,1:sample_len);
% 
% figure('Position',[50 300 1800 1000])
% for k = 1:3
%     for t = 0:length(1:ROC_dt:t_total-2*ROC_dt)
%         if k == 1
%             r1_vect = reshape(hist_vect(:,t+1,(num_trials/4 -sample_len+ 1):num_trials/4),[length(N-num_VTA+1:N) sample_len]);
%             title_1 = ['AUC for first ' num2str(sample_len) ' trials (before learning)'];
%         elseif k == 2
%             r1_vect = reshape(hist_vect(:,t+1,(3*num_trials/4 -sample_len+ 1):3*num_trials/4 ),[length(N-num_VTA+1:N) sample_len]);
%             title_1 = ['AUC for ' num2str(sample_len) ' trials after first learning (blocking)'];
%         elseif k == 3
%             r1_vect = reshape(hist_vect(:,t+1,(num_trials-sample_len + 1):num_trials),[length(N-num_VTA+1:N) sample_len]);
%             title_1 = ['AUC for ' num2str(sample_len) ' trials after learning via unblocking'];
%         end
%         for i = 1:length(N-num_VTA+1:N)
%             b_line = baseline(i,:)';
%             r1_v = r1_vect(i,:)';
%             bins = min([b_line r1_v]):max([b_line r1_v])+1;
%             prob_pos = zeros(1,length(bins));
%             prob_neg = zeros(1,length(bins));
%             for j = 1:length(bins)
%                 thresh = bins(j);
%                 prob_pos(j) = sum(r1_v>=thresh)/length(r1_v);
%                 prob_neg(j) = sum(b_line>=thresh)/length(b_line);
%             end
%             [prob_neg_sorted,idx] = sort(prob_neg);
%             prob_pos_sorted = prob_pos(idx);
%             
%             auc_plot(i,t+1) = trapz(prob_neg_sorted, prob_pos_sorted);
%         end
%     end
%     subplot(1,3,k)
%     imagesc(auc_plot)
%     caxis([0 1])
%     title(title_1)
%     xlabel('Time bin # (50ms time bins)')
%     ylabel('Neuron #')
%     mycolormap = customcolormap([0 0.425 0.5 0.575 1], [1 1 0;0 0 0;0 0 0;0 0 0;0 1 1]);
%     colorbar;
%     colormap(mycolormap)
% end