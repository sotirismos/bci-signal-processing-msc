% Script Demonstrating the SP_based Feature enginnering step AND the subsequent DataLearning procedure  
% in the case of Multicannel iEEG data  (76 sensor x 4000 samples) 
% and with the scope(s) of automated state-detection (two-class problem)  PREICTAL vs ICTAL
%                AND identifying the focus of epileptic-seizures  
%            available are 8 exemplars from each state 


clear all, close all
X=load('sz1_ict.dat');Xpre=load('sz1_pre.dat');ICTAL_DATA(:,:,1)=X;PRE_DATA(:,:,1)=Xpre;
X=load('sz2_ict.dat');Xpre=load('sz2_pre.dat');ICTAL_DATA(:,:,2)=X;PRE_DATA(:,:,2)=Xpre;
X=load('sz3_ict.dat');Xpre=load('sz3_pre.dat');ICTAL_DATA(:,:,3)=X;PRE_DATA(:,:,3)=Xpre;
X=load('sz4_ict.dat');Xpre=load('sz4_pre.dat');ICTAL_DATA(:,:,4)=X;PRE_DATA(:,:,4)=Xpre;
X=load('sz5_ict.dat');Xpre=load('sz5_pre.dat');ICTAL_DATA(:,:,5)=X;PRE_DATA(:,:,5)=Xpre;
X=load('sz6_ict.dat');Xpre=load('sz6_pre.dat');ICTAL_DATA(:,:,6)=X;PRE_DATA(:,:,6)=Xpre;
X=load('sz7_ict.dat');Xpre=load('sz7_pre.dat');ICTAL_DATA(:,:,7)=X;PRE_DATA(:,:,7)=Xpre;
X=load('sz8_ict.dat');Xpre=load('sz8_pre.dat');ICTAL_DATA(:,:,8)=X;PRE_DATA(:,:,8)=Xpre;
Fs=400, [Nsensors,Ntime]=size(X) %76 sensors x 4000 samples

figure(1),clf,subplot(1,3,1),imagesc(Xpre),xlabel('time'),ylabel('sensor no'),title('pre-ictal activity')
subplot(1,3,2),imagesc(X),xlabel('time'),ylabel('sensor no'),title('ictal activity')

% compute for each electrode the baseline mean & std of brain activity
Baseline_means=mean(Xpre,2);Baseline_stds=std(Xpre,[],2);
% use them for estimating zscore during ictal-period
subplot(1,3,3),imagesc( (X-repmat(Baseline_means,1,Ntime))./repmat(Baseline_stds,1,Ntime)), % z(t)=(xi(t)-mean(xpre(t))) /std(xpre(t))
xlabel('time'),ylabel('sensor no'),title('Rel.Increase with respect to pre-ictal'),colorbar


%figure,for i_event=1:8,subplot(2,4,i_event),imagesc(PRE_DATA(:,:,i_event)), caxis([-max(abs(PRE_DATA(:)))  max(abs(PRE_DATA(:)))]), end  
%figure,for i_event=1:8,subplot(2,4,i_event),imagesc(ICTAL_DATA(:,:,i_event)), caxis([-max(abs(ICTAL_DATA(:)))  max(abs(ICTAL_DATA(:)))]), end  
figure(2),clf
for i_event=1:8, Xpre=PRE_DATA(:,:,i_event); X=ICTAL_DATA(:,:,i_event);
Baseline_means=mean(Xpre,2);Baseline_stds=std(Xpre,[],2);
subplot(2,4,i_event),imagesc( (X-repmat(Baseline_means,1,Ntime))./repmat(Baseline_stds,1,Ntime)),title(strcat('event-',num2str(i_event))),colorbar,end


%%% PART-II
%% computing Hjorth descriptors in both conditions ('pre-ictal', 'ictal') for every sensor and for all 8 events
%Activity=var(X,[],2); % Activity is defined as the variance of each trace
%DIFF1=diff(X,1,2); % 1st derivative along the 2nd (time-dimension)
%DIFF2=diff(X,2,2); % 2nd derivative along the 2nd (time-dimension)
%Mobility=std(DIFF1,[],2)./std(X,[],2); % Mobility is the ratio of (std of fist temporal derivative-signal) to ( std of the signal)
%Complexity= (std(DIFF2,[],2)./std(DIFF1,[],2))./Mobility; % Complexity is the ratio of (Mobility of 1st derivative-signal) to (Mobility of signal)
  
H_FEATURES_PRE=[]; % 76 sensors x 3 features x 8 events
for i_event=1:8, Xpre=PRE_DATA(:,:,i_event); 
Sensor_Activity=var(Xpre,[],2);DIFF1=diff(Xpre,1,2);DIFF2=diff(Xpre,2,2);
Sensor_Mobility=std(DIFF1,[],2)./std(Xpre,[],2); 
Sensor_Complexity= (std(DIFF2,[],2)./std(DIFF1,[],2))./Sensor_Mobility;
H_FEATURES_PRE(:,:,i_event)= [Sensor_Activity Sensor_Mobility Sensor_Complexity]; end


H_FEATURES_ICTAL=[];for i_event=1:8, X=ICTAL_DATA(:,:,i_event); 
Sensor_Activity=var(X,[],2);DIFF1=diff(X,1,2);DIFF2=diff(X,2,2);
Sensor_Mobility=std(DIFF1,[],2)./std(X,[],2); 
Sensor_Complexity= (std(DIFF2,[],2)./std(DIFF1,[],2))./Sensor_Mobility;
H_FEATURES_ICTAL(:,:,i_event)= [Sensor_Activity Sensor_Mobility Sensor_Complexity]; end



AFVs=[];for ii=1:76, AFVs=[AFVs,squeeze(H_FEATURES_ICTAL(ii,:,:))'];end % ICTAL state
BFVs=[];for ii=1:76, BFVs=[BFVs,squeeze(H_FEATURES_PRE(ii,:,:))'];end  % PRE-ICTAL state
GROUP=[ones(1,8),zeros(1,8)]; % class-labels
XX=[AFVs;BFVs]'; % the whole-set of 16 Feature-Vectors (8 per state) with each FV containing 228 attributes (=3attributes per sensor), i.e. 228=3x76 

% Feature-Ranking for identifying the most infurmative features: highly-'discriminative' sensors/attributes  
[RANK,Z]=rankfeatures(XX,GROUP,'criterion','ttest'); % RANK-list contains the features from the most informative to the least informative
%[RANK,Z]=rankfeatures(XX,GROUP,'criterion','wilcoxon'); % RANK-list contains the features from the most informative to the least informative
figure,subplot(1,2,1),plot(Z),xlabel('feature no'),ylabel('Score'),title('Feature Ordering')
subplot(1,2,2),imagesc(reshape(Z,3,76)),ylabel('attribute no'),xlabel('sensor no'),title('Feature Ordering'),colorbar

                                     % alternative method for ranking-features [RANK,Z] = randfeatures(XX,GROUP,'SubsetSize',5,'Classifier','da')
best_sensor=ceil(RANK(1)/3); 
best_attribute=[(rem(RANK(1),3)==0)*3 ]+ [rem(RANK(1),3)~=0]*[rem(RANK(1),3)];

display(strcat('optimal sensor: ',num2str(best_sensor)))
display(strcat('discriminative attribute: ',num2str(  best_attribute)))

figure,
for ii=1:8, subplot(2,4,ii),plot([1:4000],ICTAL_DATA(10,:,ii),'r',[1:4000],PRE_DATA(best_sensor,:,ii),'k'),
    ylabel(strcat('sensor:',num2str(best_sensor))),title(strcat('event-',num2str(ii))),legend('ictal')
end


% just testing the classification-performance using LDA based on the best features
C = classify(XX(RANK(1:5),:)',XX(RANK(1:5),:)',GROUP);
cp = classperf(GROUP,C); display('classification-performance:'), cp.CorrectRate

