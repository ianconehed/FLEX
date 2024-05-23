function plot_func(T,delta,l,num_columns,t_j,t_total,dt,npp,N,num_VTA,R_it,R_kt,T_pt,T_dt,W_ji,Title,t_stim,dopa,dopa_plot)
%info
num_VTA = npp;
new_rate = zeros(2*num_columns,numel(t_j));
new_rate_ele = 0;
new_rate2 = zeros(2*num_columns,numel(t_j));
new_rate_ele2 = 0;
pop_vect = 0:npp:N;
set(0,'DefaultAxesColorOrder',brewermap(12,'Paired')) 
for i = 1:2*num_columns
    for t = 1:((t_total)/dt)
        for k = pop_vect(i)+1:pop_vect(i+1)
            new_rate_ele = new_rate_ele + R_it(k,t);
            new_rate_ele2 = new_rate_ele2 + R_kt(k,t);
        end

        new_rate(i,t) = new_rate_ele/(npp);
        new_rate_ele = 0;
        new_rate2(i,t) = new_rate_ele2/(npp);
        new_rate_ele2 = 0;
    end
end


Plot_array = cell(4*num_columns,1);
for i = 1:(4*num_columns)
    if mod(i,2) == 0
        Plot_array{i,1} = new_rate(i/2,:);
    else
        Plot_array{i,1} = t_j;
    end
end

Plot_array2 = cell(4*num_columns,1);
for i = 1:(4*num_columns)
    if mod(i,2) == 0
        Plot_array2{i,1} = new_rate2(i/2,:);
    else
        Plot_array2{i,1} = t_j;
    end
end

if l == 1
    f = figure('rend','painters','pos',[100 100 1500 1000]);
    p = uipanel('Parent',f,'BorderType','none');
    p.BackgroundColor = [1 1 1];
    subplot(3,1,1,'Parent',p);
    plot(Plot_array{:,1},Plot_array2{:,1},'linewidth',4);
%     plot(Plot_array{:,1},'linewidth',4);
    axis([0.01 t_total 1 140]);
    title(Title);
    ylabel('Firing Rate (Hz)');
    xticks(t_stim-ones(1,length(t_stim)));
 


    subplot(3,1,2,'Parent',p);
    plot(t_j,T_dt, 'b--',t_j,T_pt, 'r',t_j,(((T_pt.*(dopa/num_VTA)-T_dt.*(dopa_plot/num_VTA)))*3)+70,'g--','linewidth',4);
    axis([0.01 t_total 1 140]);
    xticks(t_stim-ones(1,length(t_stim)));
    xlabel('Time(ms)');
    ylabel('Eligibility Trace (AU)');
    
    subplot(3,1,3,'Parent',p);
    plot(1:dt:t_total,mean(R_it(N-num_VTA+1:N,:)),'linewidth',4);
    hold on
    plot(1:dt:t_total,mean(R_kt(npp+1:2*npp,:)),'linewidth',4);
    hold on
    plot(1:dt:t_total,mean(R_kt(N-num_VTA+1:N,:)),'linewidth',4);
    hold off
    axis([0.01 t_total 1 25]);
    xticks(t_stim-ones(1,length(t_stim)));
    xlabel('Time(ms)');
    ylabel('Eligibility Trace (AU)');
    
%     sp4 = subplot(1,2,2,'Parent',p);
%     hold on
%     for i = (N-num_VTA+1):N
%         plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),running_mean_R_it(i,:)*1000)
%     end
%     view(0,85)
%     hold off
%     xlim([0 t_total])
%     xlabel('Time(ms)')
%     ylabel('Neuron number')
%     zlabel('Firing rate (Hz)')
%     title(['VTA neurons at trial ' num2str(l) ' averaged for ' num2str(25*(l>=25) + l*(l<25)) ' trials'])
%     
%     set(sp4, 'DefaultAxesColorOrder', colororder(white(100)-1))


    p.FontSize = 34;
else
    subplot(3,1,1);
    plot(Plot_array{:,1},Plot_array2{:,1},'linewidth',4);
%     plot(Plot_array{:,1},'linewidth',4);
    axis([0.01 t_total 1 140]);
    title(Title);
    ylabel('Firing Rate (Hz)');
    xticks(t_stim-ones(1,length(t_stim)));
 


    subplot(3,1,2);
    plot(t_j,T_dt, 'b--',t_j,T_pt, 'r',t_j,(((T_pt.*(dopa/num_VTA)-T_dt.*(dopa_plot/num_VTA)))*3)+70,'g--','linewidth',4);
    axis([0.01 t_total 1 140]);
    xticks(t_stim-ones(1,length(t_stim)));
    xlabel('Time(ms)');
    ylabel('Eligibility Trace (AU)');
    int_sum = sum((((T_pt.*(dopa/num_VTA)-T_dt.*(dopa_plot/num_VTA)))/3));
    text(100,120,num2str(int_sum))
    
    subplot(3,1,3);
    plot(1:dt:t_total,mean(R_it(N-num_VTA+1:N,:)),'linewidth',4);
    hold on
    plot(1:dt:t_total,mean(R_kt(npp+1:2*npp,:)),'linewidth',4);
    hold on
    plot(1:dt:t_total,mean(R_kt(N-num_VTA+1:N,:)),'linewidth',4);
    hold off
    axis([0.01 t_total 1 25]);
    xticks(t_stim-ones(1,length(t_stim)));
    xlabel('Time(ms)');
    ylabel('Eligibility Trace (AU)');
    
%     sp4 = subplot(1,2,2);
%     hold on
%     for i = (N-num_VTA+1):N
%         plot3(1:dt:t_total,i*ones(size(1:dt:t_total)),running_mean_R_it(i,:))
%     end
%     view(0,85)
%     hold off
%     xlim([0 t_total])
%     xlabel('Time(ms)')
%     ylabel('Neuron number')
%     zlabel('Firing rate (Hz)')
%     title(['VTA neurons at trial ' num2str(l) ' averaged for ' num2str(25*(l>=25) + l*(l<25)) ' trials'])
%     
%     set(sp4, 'DefaultAxesColorOrder', colororder(white(100)-1))
end
drawnow
end


