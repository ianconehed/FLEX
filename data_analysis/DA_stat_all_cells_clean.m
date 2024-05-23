%get statistics from Amo data
%First time learining data from Amo - assumes data in same directory as .m
%file. Otherwise need to set path for the .mat file
%Only first 40 trials, or less if there are less type 1 trials are used
%This is because ther seems to be a habituation of US response for a larger
%number of trials

%Cell_number=440; %437-440,444-446 (7 animals) 

Cell_number_list=[437:440,444:446]; % make sure these data files are in this directory

Num_days=length(Cell_number_list);
number_of_days=10;
mean_dyn=zeros(number_of_days,9001);
cs_peak=zeros(number_of_days,Num_days);
us_peak=zeros(number_of_days,Num_days);
DA_integral=zeros(number_of_days,Num_days);



for dd=1:Num_days
    Cell_number=Cell_number_list(dd)
    for ii=1:number_of_days
        mat_file=['FirstTimeLearning_',num2str(Cell_number),'_day',num2str(ii),'.mat']; %this sets and loads Amo .mat files
        load(mat_file);
        vn1=['Trial_number_lick_',num2str(ii),'=Trial_number_lick;'];
        eval(vn1);
        vn2=['DeltaF_licktrial_',num2str(ii),'=DeltaF_licktrial;'];
        eval(vn2);
        kk=min(Trial_number_lick(1),40); %uses 40 or the number of type 1 trials, if they are smaller than 40
        mean_dyn(ii,:)=mean(DeltaF_licktrial(1:kk,:));
        cs_peak(ii,dd)=max(mean_dyn(ii,2000:4000));
        us_peak(ii,dd)=max(mean_dyn(ii,5500:9000));
        DA_integral(ii,dd)=sum(mean_dyn(ii,2000:end));
        %figure(ii);
        %imagesc(DeltaF_licktrial(1:kk,:));
        
    end
    
    figure(30);
    
    plot(DA_integral(:,dd),'+','LineWidth',2);
   
    ylabel('DA integral [AU]','FontSize',18)
    xlabel('Day number','FontSize',18)
    hold on
    eval(['DA_integral_',num2str(Cell_number),'=DA_integral;'])

    eval(['us_peak_',num2str(Cell_number),'=us_peak;'])
    eval(['cs_peak_',num2str(Cell_number),'=cs_peak;'])
    eval(['DA_integral_',num2str(Cell_number),'=DA_integral;'])
    
    eval(['save ','processed_data_',num2str(Cell_number)])
end

figure(30)
hold on
DA_integral_mean=mean(DA_integral');
plot(DA_integral_mean,'LineWidth',2);


Day_1_2_data=[DA_integral(1,:),DA_integral(2,:)];
Day_3_4_data=[DA_integral(3,:),DA_integral(4,:)];
Day_8_10_data=[DA_integral(8,:),DA_integral(9,:),DA_integral(10,:)];
Day_8_9_data=[DA_integral(8,:),DA_integral(9,:)];
DA_1_2_mean=mean(Day_1_2_data);
DA_1_2_stderr=std(Day_1_2_data)/sqrt(length(Day_1_2_data));
DA_1_2_mean_low=DA_1_2_mean-DA_1_2_stderr;
DA_1_2_mean_high=DA_1_2_mean+DA_1_2_stderr;

DA_3_4_mean=mean(Day_3_4_data);
DA_3_4_stderr=std(Day_3_4_data)/sqrt(length(Day_3_4_data));
DA_3_4_mean_low=DA_3_4_mean-DA_3_4_stderr;
DA_3_4_mean_high=DA_3_4_mean+DA_3_4_stderr;

DA_8_10_mean=mean(Day_8_10_data);
DA_8_10_stderr=std(Day_8_10_data)/sqrt(length(Day_8_10_data));
DA_8_10_mean_low=DA_8_10_mean-DA_8_10_stderr;
DA_8_10_mean_high=DA_8_10_mean+DA_8_10_stderr;

Stat_mean=[DA_1_2_mean,DA_3_4_mean,DA_8_10_mean];
Stat_low=[DA_1_2_mean_low,DA_3_4_mean_low,DA_8_10_mean_low];
Stat_high=[DA_1_2_mean_high,DA_3_4_mean_high,DA_8_10_mean_high];
figure(50)
x=1:3;
bar(x,Stat_mean)
hold on
er= errorbar(x,Stat_mean,Stat_low,Stat_high);
er.LineWidth = 3; 
er.Color=[0 0 0]

[p_1_2v3_4,h_1_2v3_4]=ranksum(Day_1_2_data,Day_3_4_data);

if h_1_2v3_4==1
    ['RankSum statistics of days 1,2 vs. days 3,4 is p=',num2str(p_1_2v3_4)]
else
    ['RankSum statistics of days 1,2 vs. days 3,4 is not significant']
end

[p_1_2v8_10,h_1_2v8_10]=ranksum(Day_1_2_data,Day_8_10_data);

if h_1_2v8_10==1
    ['RankSum statistics of days 1,2 vs. days 8-10 is p=',num2str(p_1_2v8_10)]
else
    ['RankSum statistics of days 1,2 vs. days 8-10 is not significant']
end
    
%Also works for days 8,9 only. However not enough data for one day at a
%time

[p_1_2v8_10,h_1_2v8_9]=ranksum(Day_1_2_data,Day_8_9_data);

%Also works with paired ttest for 1,2 vs 3,4 
[pt,ht]=ttest(Day_1_2_data,Day_3_4_data);
%-not for 1,2 vs 8,9 - which is less important as it should theoretically
% come down at some point anyway.
% No real way to justify normality of data though, so ranksum is preferable 



%Another way of doing the stats - average over the days, this reduces
%variablity but also reduces the number of data points
Da_1_2_mean=mean(DA_integral(1:2,:));
Da_3_4_mean=mean(DA_integral(3:4,:));
[pm,hm]=ranksum(Da_1_2_mean,Da_3_4_mean);
%This is still signficant
Da_8_10_mean=mean(DA_integral(8:10,:));
[pm1,hm1]=ranksum(Da_1_2_mean,Da_8_10_mean);
% This is not though, so again, it comes down
[pm2,hm2]=ranksum(Da_3_4_mean,Da_8_10_mean);
% but not significant


