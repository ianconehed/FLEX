%% average weights plotting
figure('Position',[1200 300 500 500])
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

lengthie = 30;
minus_trial = ceil(num_trials/3); 


mean_VTA_plot = mean(R_it(N-num_VTA+1:N,:));
mean_VTA_plot_minus = mean(R_it_l(N-num_VTA+1:N,:,minus_trial),1);
% mean_R_it_after = mean(R_it_l(:,:,(num_trials-lengthie + 1):num_trials),3);
mean_R_it_after = mean(R_it_l(:,:,num_trials-lengthie+1:num_trials),3);
mean_R_it_after_learning = mean(R_it_l(:,:,3*num_trials/4-lengthie+1:3*num_trials/4),3);
mean_R_it_before = mean(R_it_l(:,:,1:lengthie),3);
% mean_R_it_before = mean(R_it_l(:,:,2*num_trials/4 - lengthie + 1:2*num_trials/4),3);
mean_mean_R_it_after = mean(mean_R_it_after(N-num_VTA+1:N,:));
mean_mean_R_it_after_learning = mean(mean_R_it_after_learning(N-num_VTA+1:N,:));
mean_mean_R_it_before = mean(mean_R_it_before(N-num_VTA+1:N,:));
mean_R_it = mean(R_it_l,3);

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

%% spike raster

figure('Position',[300 300 1000 1000])
subplot(2,1,1)
LineFormat = struct();
LineFormat.Color = 'yellow';
plotSpikeRaster(sc_R_it(:,:,num_trials) == 1,'LineFormat',LineFormat);
set(gca,'color',[0 0 0])
ylabel('Neuron number')
xlabel('Time(ms)')
title('spike raster  (200-300 are VTA neurons)')
subplot(2,1,2)
plot(1:dt:t_total,mean_VTA_plot)
ylabel('Firing rate')
xlabel('Time(ms)')
title('Mean VTA neuron firing rate, last trial')
xlim([0.01 t_total]);
ylim([0.01 max(mean_VTA_plot) + 5]);

figure('Position',[300 300 1000 1000])
subplot(2,1,1)
LineFormat = struct();
LineFormat.Color = 'yellow';
plotSpikeRaster(sc_R_it(:,:,minus_trial) == 1,'LineFormat',LineFormat);
set(gca,'color',[0 0 0])
ylabel('Neuron number')
xlabel('Time(ms)')
title('spike raster  (200-300 are VTA neurons)')
subplot(2,1,2)
plot(1:dt:t_total,mean_VTA_plot_minus)
ylabel('Firing rate')
xlabel('Time(ms)')
title('Mean VTA neuron firing rate, last trial')
xlim([0.01 t_total]);
ylim([0.01 max(mean_VTA_plot_minus) + 5]);


%% individual averaged plotting

figure('Position',[50 300 1800 1000])
subplot(1,3,1)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_before(i,:))
end
view(0,85)
xlim([0 t_total])
xlabel('Time(ms)')
ylabel('Neuron number')
zlabel('Firing rate (Hz)')
title(['VTA neurons, averaged for ' num2str(lengthie) ' trials before learning'])

colororder(white(100)-1)

subplot(1,3,2)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_after_learning(i,:))
end
view(0,85)
xlim([0 t_total])
xlabel('Time(ms)')
ylabel('Neuron number')
zlabel('Firing rate (Hz)')
title(['VTA neurons, averaged for ' num2str(lengthie) ' trials after learning (blocking)'])

colororder(white(100)-1)

subplot(1,3,3)
hold on
for i = (N-num_VTA+1):N
    plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),mean_R_it_after(i,:))
end
view(0,85)
xlim([0 t_total])
xlabel('Time(ms)')
ylabel('Neuron number')
zlabel('Firing rate (Hz)')
title(['VTA neurons, averaged for ' num2str(lengthie) ' trials after unblocking'])

colororder(white(100)-1)


%% mean plotting

figure
subplot(3,1,1)
plot(1:dt:t_total,mean_mean_R_it_before)
title(['mean over all neurons and'  num2str(lengthie) 'trials before learning'])
xlabel('time')
ylabel('Firing rate')
subplot(3,1,2)
plot(1:dt:t_total,mean_mean_R_it_after_learning)
title(['mean over all neurons and last' num2str(lengthie) 'trials of learning'])
xlabel('time')
ylabel('Firing rate')
subplot(3,1,3)
plot(1:dt:t_total,mean_mean_R_it_after)
title(['mean over all neurons and' num2str(lengthie) 'trials of no CS'])
xlabel('time')
ylabel('Firing rate')




%% AUC plotting


baseline = hist_vect1(:,1:lengthie);

figure('Position',[50 300 1800 1000])
for k = 1:3
    for t = 0:length(1:ROC_dt:t_total-2*ROC_dt)
        if k == 1
            r1_vect = reshape(hist_vect(:,t+1,(num_trials/4 -lengthie+ 1):num_trials/4),[length(N-num_VTA+1:N) lengthie]);
            title_1 = ['AUC for first ' num2str(lengthie) ' trials (before learning)'];
        elseif k == 2
            r1_vect = reshape(hist_vect(:,t+1,(3*num_trials/4 -lengthie+ 1):3*num_trials/4 ),[length(N-num_VTA+1:N) lengthie]);
            title_1 = ['AUC for ' num2str(lengthie) ' trials after first learning (blocking)'];
        elseif k == 3
            r1_vect = reshape(hist_vect(:,t+1,(num_trials-lengthie + 1):num_trials),[length(N-num_VTA+1:N) lengthie]);
            title_1 = ['AUC for ' num2str(lengthie) ' trials after learning via unblocking'];
        end
        for i = 1:length(N-num_VTA+1:N)
            b_line = baseline(i,:)';
            r1_v = r1_vect(i,:)';
            bins = min([b_line r1_v]):max([b_line r1_v])+1;
            prob_pos = zeros(1,length(bins));
            prob_neg = zeros(1,length(bins));
            for j = 1:length(bins)
                thresh = bins(j);
                prob_pos(j) = sum(r1_v>=thresh)/length(r1_v);
                prob_neg(j) = sum(b_line>=thresh)/length(b_line);
            end
            [prob_neg_sorted,idx] = sort(prob_neg);
            prob_pos_sorted = prob_pos(idx);
            
            auc_plot(i,t+1) = trapz(prob_neg_sorted, prob_pos_sorted);
        end
    end
    subplot(1,3,k)
    imagesc(auc_plot)
    caxis([0 1])
    title(title_1)
    xlabel('Time bin # (50ms time bins)')
    ylabel('Neuron #')
    mycolormap = customcolormap([0 0.425 0.5 0.575 1], [1 1 0;0 0 0;0 0 0;0 0 0;0 1 1]);
    colorbar;
    colormap(mycolormap)
end