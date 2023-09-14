%get statistics from Amo data
%First time learining
%Only first 40 trials, or less if there are less type 1 trials are used
%This is because ther seems to be a habituation of US response for a larger
%number of trials
% single animal comparison between days, based on multiple trials in each
% day

Cell_number=446; %437-440,444-446 (7 animals) 




number_of_days=10;
mean_dyn=zeros(number_of_days,9001);
cs_peak=zeros(number_of_days,1);
us_peak=zeros(number_of_days,1);
DA_integral=zeros(number_of_days,1);

for ii=1:number_of_days
    mat_file=['FirstTimeLearning_',num2str(Cell_number),'_day',num2str(ii),'.mat'];
    load(mat_file);
    vn1=['Trial_number_lick_',num2str(ii),'=Trial_number_lick;'];
    eval(vn1);
    vn2=['DeltaF_licktrial_',num2str(ii),'=DeltaF_licktrial;'];
    eval(vn2);
    kk=min(Trial_number_lick(1),40); %uses 40 or the number of type 1 trials, if they are smaller than 40
    mean_dyn(ii,:)=mean(DeltaF_licktrial(1:kk,:));
    eval(['Dyn_day_trial_',num2str(ii),'=DeltaF_licktrial(1:',num2str(kk),',:);'])
    eval(['Int_day_trail_',num2str(ii),'=mean(transpose(Dyn_day_trial_',num2str(ii),'));'])
    DA_integral(ii)=sum(mean_dyn(ii,2000:end));
    %figure(ii);
    %imagesc(DeltaF_licktrial(1:kk,:));
    
end 
figure(20)
plot(mean_dyn','LineWidth',2);
legend(num2str([1:number_of_days]'));
eval(['mean_dyn_',num2str(Cell_number),'=mean_dyn;'])
eval(['mean_dyn_',num2str(Cell_number),'=mean_dyn;'])
figure(30);
%
plot(DA_integral,'+','LineWidth',2);
title('DA integral');
% hold on
[a13,b13]=ranksum(Int_day_trail_1,Int_day_trail_3)
sum(Int_day_trail_1)<sum(Int_day_trail_3)
[a14,b14]=ranksum(Int_day_trail_1,Int_day_trail_4)
sum(Int_day_trail_1)<sum(Int_day_trail_4)
[a23,b23]=ranksum(Int_day_trail_2,Int_day_trail_3)
sum(Int_day_trail_2)<sum(Int_day_trail_3)
[a24,b24]=ranksum(Int_day_trail_2,Int_day_trail_4)
sum(Int_day_trail_2)<sum(Int_day_trail_4)
[a18,b18]=ranksum(Int_day_trail_1,Int_day_trail_8)
sum(Int_day_trail_1)<sum(Int_day_trail_8)
[a19,b19]=ranksum(Int_day_trail_1,Int_day_trail_9)
sum(Int_day_trail_1)<sum(Int_day_trail_9)
[a28,b28]=ranksum(Int_day_trail_2,Int_day_trail_8)
sum(Int_day_trail_2)<sum(Int_day_trail_8)
[a29,b29]=ranksum(Int_day_trail_2,Int_day_trail_9)
sum(Int_day_trail_2)<sum(Int_day_trail_9)




% eval(['DA_integral_',num2str(Cell_number),'=DA_integral;'])
% 
% eval(['us_peak_',num2str(Cell_number),'=us_peak;'])
% eval(['cs_peak_',num2str(Cell_number),'=cs_peak;'])
% eval(['DA_integral_',num2str(Cell_number),'=DA_integral;'])
% 
% eval(['save ','processed_data_',num2str(Cell_number)])
