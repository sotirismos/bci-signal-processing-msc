% script for demonstrating the use Vector-Ordering for detecting outlying Single-Trial Respones.
% A set of 110 Single-Trial Visual Evoked Responses, contaminated by ''artifacts''(i.e. outliers), is used.
% The task is to remove the outliers and study their impact in terms of SNR (signal-to-noise ratio).

clear all, close all
load VEPS
% Get the first 50 single-trial responses out of the 110.
%veps=veps(1:50,:);
[Ntrials,Nsamples]=size(veps);

% sets the sampling frequency to 1000 Hz, meaning each sample represents 1/1000 of a second.
Fs=1000;
% creates a vector t from 1 to the number of samples, which serves as the x-axis in subsequent plots.
t=1:Nsamples;

% plot(t,mean(veps)) computes and plots the average response across all trials as a function of sample number.
% The x-axis is labeled “sample no”, the title is set to “averaged response”, and the y-axis is labeled with its units (arbitrary units in volts).
% figure,subplot(2,1,1),plot(t, mean(veps)), xlabel('sample no'), title('averaged response'), ylabel('a.u.(volts)')

% here, the x-axis is converted into seconds by multiplying t by the sample period (1/Fs).
% subplot(2,1,2), plot(t*(1/Fs),mean(veps)), xlabel('time (s)')

% the detrend function removes any constant (DC) offset (or even linear trends) from the data.
VEPS=detrend(veps')';

% plot(t,veps) plots all original (non-detrended) single-trial waveforms over time.
% figure, subplot(1,2,1), plot(t,veps), hold, plot(t,mean(veps),'k','linewidth',2), hold
% xlabel('sample no'), title('single-trial responses - original level'), ylabel('a.u.(volts)')

% plots the detrended waveforms. This side-by-side comparison shows how the baseline has been removed (“dc-offset” removed).
% subplot(1,2,2),plot(t,VEPS),hold,plot(t,mean(VEPS),'k','linewidth',2),hold
% xlabel('sample no'),title('after applying dc-offset'),ylabel('a.u.(volts)')


% measuring the SNR (after removing the DC-level)
% call a custom function snr_sample (snr_sample.m) to compute the estimated signal power (sp) and noise power (np)
[sp,np] = snr_sample(VEPS); 

% computes the SNR for an individual trial as the ratio of signal to noise power.
trial_SNR = sp/np; 
% Calculates the SNR for the averaged response by multiplying the single-trial SNR by the number of trials.
ave_SNR = Ntrials*trial_SNR;

disp('signal-to-noise ratio per trial:'), trial_SNR;
disp('SNR of the averaged response:'), ave_SNR;

% using 'pairwise-distances' to reveal the presence of outliers

% squareform converts the condensed distance vector into a full symmetric matrix (Dmatrix), where each element (i, j) represents the distance between trial i and trial j.
Dmatrix = squareform(pdist(VEPS));
disp('size of distance matrix:'), size(Dmatrix)

% figure, subplot(1,2,1), plot(t,VEPS),title('single-trial waveforms'), xlabel('time (sample no)')
% subplot(2,2,2), imagesc(Dmatrix),axis square, title('Distance-Matrix'), xlabel('trial no'), ylabel('trial no')

% compute an aggregate distance for each trial by summing the distances of that trial to all other trials.
Dist_Score = sum(Dmatrix); 
% subplot(2,2,4), stem(Dist_Score,'o-k','markerfacecolor','r'), grid, xlim([1 Ntrials]) 
% xlabel('trial no'), ylabel('aggregate distance'), hold


% ranking the waveforms based on aggregate-distance

% The sort function sorts the aggregate distance scores in ascending order.
[Ranked_Dist_Score, sorted_list] = sort(Dist_Score);

% Plots the sorted (ranked) aggregate distance scores
figure,subplot(1,2,1),plot(Ranked_Dist_Score,'k-o','markerfacecolor','r'),

disp('the candidate outliers:')
% Selects the last six (or three) entries from sorted_list in descending order
outlier_list = num2cell(sorted_list(end:-1:end-5));

text((Ntrials:-1:Ntrials-5)+1,Ranked_Dist_Score(end:-1:end-5),outlier_list)
%subplot(1,2,2),imagesc(VEPS(sorted_list,:)),xlabel('time (sample no)'),ylabel('ranked-waveforms'),title('sorted ST-waveforms')


% removing the outliers

outlier_list=sorted_list(end:-1:end-5);

% kept_list holds the indices of trials that are considered “clean.”
kept_list=setdiff(1:Ntrials,outlier_list);

% redefines (or confirms) the outlier_list as the set of six trials with the highest aggregate distance scores.
%figure, subplot(1,2,1), plot(1:Nsamples,VEPS(kept_list,:),'k',1:Nsamples,VEPS(outlier_list,:),'r'),
%title('single-trial waveforms'), xlabel('time (sample no)')

%subplot(1,2,2), plot(1:Nsamples,mean(VEPS(kept_list,:)),'b',1:Nsamples,mean(VEPS),'r'),
%legend('selective-averaging','ensemble-averaging'),xlabel('time (sample no)')

% the script calls the snr_sample function again, this time using only the kept trials (after removing the outliers) to calculate the new signal and noise power.
[sp,np] = snr_sample(VEPS(kept_list,:)); 
selective_trial_SNR=sp/np; 
disp('relative increase in trial-level SNR:')
(selective_trial_SNR-trial_SNR)/trial_SNR
              

