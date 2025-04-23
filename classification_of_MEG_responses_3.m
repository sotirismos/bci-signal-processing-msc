% The data are coming from an Auditory stimulation paradigm. 
% In two distinct files, the signals corresponding to right/left auditory
% cortex are included.
% Additional 'control' data (with no stimulus) have also beed included
% and utilized for comparison purposes.

% NOTE: it is assumed that the directory : Auditory_MEG_responses has been added to the path.
% Otherwise uncomment the following line and the line at the end of this script
% addpath stimulation, addpath control

clear all, close all

%% Step-1
% We load 2 MEG long-signals (1 for each hemisphere) and 1 sequence containing the times when the stimulus was delivered
load stim_times , load subj3_left, load subj3_right
Fs=625; % sampling frequency

% quick data sanity check
fprintf('Sampling frequency: %d Hz\n', Fs);
fprintf('Continuous recording len: %d samples (%.2f s)\n', numel(subj3_left), numel(subj3_left)/Fs);
fprintf('Number of stimuli: %d\n', numel(stim_times(1:end-1)));

%% Visual inspection of stimulation recording
% The right hemisphere signal is offset by adding 5 (i.e. subj3_right+5)
% for clarity—offsetting prevents overlap and makes comparisons easier.
% Stimulus onsets are overlaid as black dots using the stim_times for the x-axis 
% and a constant negative value (-2) for the y-axis, ensuring they are
% visible in the context of the MEG signals.

figure(1),clf,subplot(2,1,1),plot(1:numel(subj3_left), subj3_left, 1:numel(subj3_right), subj3_right+5, stim_times, -2*ones(1,numel(stim_times)), 'k.'),
legend('left','right','stim onset'),xlabel('time (sample no)'),title('stimulation')
 
%% Load and plot spontaneous/control activity
load subj3_control_left, load subj3_control_right
subplot(2,1,2),plot(1:numel(subj3_control_left), subj3_control_left,1:numel(subj3_control_right),subj3_control_right+5),
legend('left','right'),xlabel('time (sample no)'),title('spontaneous activity')

%% Step-2  
% Frequency-Content Representation by means of Power-SPectral-Density (PSD)
% TODO: Analyze how pwelch works and the parameters used.
% TODO: analyze the y-axis (Power/Frequency (dB/Hz)
figure(2);
subplot(2,2,1), pwelch(subj3_left,1024,500,1024,Fs,'onesided'), title('PSD: stimulation'),legend('left')
subplot(2,2,2), pwelch(subj3_right,1024,500,1024,Fs,'onesided'), title('PSD: stimulation'),legend('right')
subplot(2,2,3), pwelch(subj3_control_left,1024,500,1024,Fs,'onesided'), title('PSD: spontaneous'),legend('left')
subplot(2,2,4), pwelch(subj3_control_right,1024,500,1024,Fs,'onesided'), title('PSD: spontaneous'),legend('right')
% alternatively we may use sptool to derive and compare the spectra (TODO: analyze sptool)           
              
%% Step-3  
% Segmenting STIMULATION data into trials.
% Using the following simple m-file signal_to_trials.m 
% we extract segments from the continuous signal (stimulation condition)
% around the times of the stimules onset (stim_times variable)
% TODO: Analyze M100
[left_trials,~]=signal_to_trials(subj3_left,stim_times,100,300); % 100 samples before the stim-onset and 300 samples after
[right_trials,t]=signal_to_trials(subj3_right,stim_times,100,300); % vector t keeps the latencies [-100 : 300] 

figure(3);
subplot(2,2,1), plot(t*(1/Fs),left_trials(1:10:end,:),'b'),title('left-hemisphere single-trials'),xlabel('time(s)')
subplot(2,2,3), plot(t*(1/Fs),mean(left_trials),'b'),title('M100 (averaged) response'),xlabel('time(s)'),ylabel('a.u.(Tesla)'),legend('left')
subplot(2,2,2), plot(t*(1/Fs),right_trials(1:10:end,:),'r'),title('right-hemisphere single-trials'),xlabel('time(s)')
subplot(2,2,4), plot(t*(1/Fs),mean(right_trials),'r'),title('M100 (averaged) response'),xlabel('time(s)'),ylabel('a.u.(Tesla)'),legend('right')

figure(4);
subplot(1,2,1), imagesc(left_trials), colormap(jet),title('left hemisphere responses'),ylabel('trial no'),xlabel('time(sample no)'),clim([-2 2])
subplot(1,2,2), imagesc(right_trials),title('right hemisphere responses'),ylabel('trial no'),xlabel('time(sample no)'),clim([-2 2])


%% Step-4
% Segmenting CONTROL data into pseudo‑trials
% Using the following simple m-file signal_to_trials.m 
% we extract segments from the continuous control signal
% around the times of the stimules onset (stim_times variable)
[ctrl_left_trials,~]=signal_to_trials(subj3_control_left,stim_times,100,300); % 100 samples before the stim-onset and 300 samples after
[ctrl_right_trials,t]=signal_to_trials(subj3_control_right,stim_times,100,300); % vector t keeps the latencies [-100 : 300] 

figure(5);
subplot(2,2,1), plot(t*(1/Fs),ctrl_left_trials(1:10:end,:),'b'),title('left-hemisphere control single-trials'),xlabel('time(s)')
subplot(2,2,3), plot(t*(1/Fs),mean(ctrl_left_trials),'b'),title('M100 (averaged) response - control'),xlabel('time(s)'),ylabel('a.u.(Tesla)'),legend('left')
subplot(2,2,2), plot(t*(1/Fs),ctrl_right_trials(1:10:end,:),'r'),title('right-hemisphere control single-trials'),xlabel('time(s)')
subplot(2,2,4), plot(t*(1/Fs),mean(ctrl_right_trials),'r'),title('M100 (averaged) response - control'),xlabel('time(s)'),ylabel('a.u.(Tesla)'),legend('right')

figure(6);
subplot(1,2,1), imagesc(ctrl_left_trials), colormap(jet),title('left hemisphere responses - control'),ylabel('trial no'),xlabel('time(sample no)'),clim([-2 2])
subplot(1,2,2), imagesc(ctrl_right_trials),title('right hemisphere responses - control'),ylabel('trial no'),xlabel('time(sample no)'),clim([-2 2])


%% Step-5  
% Band-Pass filtering of single-trials and control pseudo-trials.
% Filtering within [3-20] Hz frequency band.
% TODO: Analyze how butter works and its parameters.
[b,a]=butter(5,[3,20]/(Fs/2)) ;   % We design the filter considering the sampling frequency

%filtering left-hemisphere responses
f_left_trials=[filtfilt(b,a,left_trials')]'; % We apply the the filter (in zero-phase mode). 
figure(5);
subplot(2,2,1), plot(t*(1/Fs), f_left_trials(1,:), t*(1/Fs), left_trials(1,:),'k'),legend('filtered','original'),hold
                plot(t*(1/Fs), f_left_trials(10,:)+2, t*(1/Fs), left_trials(10,:)+2,'k'),hold,title('left hemisphere single-trials')
subplot(2,2,3), plot(t*(1/Fs), mean(f_left_trials),'b',t*(1/Fs), mean(left_trials),'k'),legend('filtered','original'),title('left hemisphere M100 (averaged) response')

%filtering right-hemisphere responses
f_right_trials=[filtfilt(b,a,right_trials')]'; % We apply the the filter (in zero-phase mode). 
subplot(2,2,2), plot(t*(1/Fs), f_right_trials(1,:), t*(1/Fs), right_trials(1,:),'k'),legend('filtered','original'),hold
                plot(t*(1/Fs), f_right_trials(10,:)+2, t*(1/Fs), right_trials(10,:)+2,'k'),hold,title('right hemisphere single-trials')
subplot(2,2,4), plot(t*(1/Fs), mean(f_right_trials),'b', t*(1/Fs), mean(right_trials),'k'),legend('filtered','original'),title('right hemisphere M100 (averaged) response')
  
figure(6);
subplot(1,2,1), imagesc(f_left_trials),colormap(jet),title('left hemisphere responses'),ylabel('trial no'),xlabel('time(sample no)'),clim([-1.7 1.7])
subplot(1,2,2), imagesc(f_right_trials),title('right hemisphere responses'),ylabel('trial no'),xlabel('time(sample no)'),clim([-1.7 1.7])

%filtering left-hemisphere control responses
ctrl_f_left_trials=[filtfilt(b,a,ctrl_left_trials')]'; % We apply the the filter (in zero-phase mode). 
figure(7);
subplot(2,2,1), plot(t*(1/Fs), ctrl_f_left_trials(1,:), t*(1/Fs), ctrl_left_trials(1,:),'k'),legend('filtered','original'),hold
                plot(t*(1/Fs), ctrl_f_left_trials(10,:)+2, t*(1/Fs), ctrl_left_trials(10,:)+2,'k'),hold,title('left hemisphere preudo single-trials (control)')
subplot(2,2,3), plot(t*(1/Fs), mean(ctrl_f_left_trials),'b',t*(1/Fs), mean(ctrl_left_trials),'k'),legend('filtered','original'),title('left hemisphere M100 (averaged) response - control')

%filtering right-hemisphere control responses
ctrl_f_right_trials=[filtfilt(b,a,ctrl_right_trials')]'; % We apply the the filter (in zero-phase mode). 
subplot(2,2,2), plot(t*(1/Fs), ctrl_f_right_trials(1,:), t*(1/Fs), ctrl_right_trials(1,:),'k'),legend('filtered','original'),hold
                plot(t*(1/Fs), ctrl_f_right_trials(10,:)+2, t*(1/Fs), ctrl_right_trials(10,:)+2,'k'),hold,title('right hemisphere pseudo single-trials (control)')
subplot(2,2,4), plot(t*(1/Fs), mean(ctrl_f_right_trials),'b', t*(1/Fs), mean(ctrl_right_trials),'k'),legend('filtered','original'),title('right hemisphere M100 (averaged) response - control')
  
  
figure(8);
subplot(1,2,1), imagesc(ctrl_f_left_trials),colormap(jet),title('left hemisphere control responses'),ylabel('trial no'),xlabel('time(sample no)'),clim([-1.7 1.7])
subplot(1,2,2), imagesc(ctrl_f_right_trials),title('right hemisphere control responses'),ylabel('trial no'),xlabel('time(sample no)'),clim([-1.7 1.7])

%% Step-6
% Feature Extraction
% Combined time-domain, spectral and wavelet features for each trial
% Outputs: featMat_stim, featMat_ctrl, labels
trials_number=length(left_trials(:,1));
td_count=4; sp_count=3; wv_count=3;
features_per_hemi=td_count + sp_count + wv_count; 
total_feats=features_per_hemi*2;
feature_matrix=zeros(trials_number*2, total_feats);
labels=[ones(trials_number,1); zeros(trials_number,1)];

% Latency masks
pre_stimulus = t<0;
post_stimulus = t>=0;
win_m100 = t>=50 & t<=140;
win_auc = t>=0  & t<=200;

row = 1;
% ----- per condition (1=stim, 2=control) -----
for cond = 1:2
    % ----- per trial (or pseudo-trial) -----
    for k=1:trials_number
        if cond==1
            x_left = f_left_trials(k,:);
            x_right = f_right_trials(k,:);
        else
            x_left = ctrl_f_left_trials(k,:);
            x_right = ctrl_f_right_trials(k,:);
        end
        hemi_data = {x_left, x_right};
        col = 1;
         % ----- per hemisphere -----
        for h=1:2
            x = hemi_data{h};

            % ----- baseline-correct -----
            x = x - mean(x(pre_stimulus));

            % ----- Time-domain -----
            % area 0–200
            feature_matrix(row,col)=trapz(x(win_auc)); col=col+1;
            % RMS
            feature_matrix(row,col)=rms(x(post_stimulus)); col=col+1;
            % Mobility, Complexity
            [~,mob,cmp] = hjorth(x(post_stimulus)); feature_matrix(row,col)=mob; col=col+1;
            feature_matrix(row,col)=cmp; col=col+1;

            % ----- Spectral -----
            bp_theta = bandpower(x,Fs,[4 7]); 
            bp_alpha = bandpower(x,Fs,[8 12]); 

            [Pxx,freq] = pwelch(x,256,128,256,Fs,'onesided');
            sc = sum(freq.*Pxx)/sum(Pxx);

            feature_matrix(row,col)=bp_theta; col=col+1;
            feature_matrix(row,col)=bp_alpha; col=col+1;
            feature_matrix(row,col)=sc; col=col+1;

            % ----- Wavelet -----
            % 3-level
            [c,l] = wavedec(x,3,'db4'); 
            d3=detcoef(c,l,3); 
            d2=detcoef(c,l,2); 
            a3=appcoef(c,l,'db4',3);
            
            feature_matrix(row,col)=sum(d3.^2); col=col+1;
            feature_matrix(row,col)=sum(d2.^2); col=col+1;
            feature_matrix(row,col)=sum(a3.^2); col=col+1;
        end
        row = row + 1;
    end
end

%% Step-7 Classification with Cross-Validation
% Prepare features and labels
X = feature_matrix;
Y = labels;

% Stratified 5-fold partition
cvp = cvpartition(Y,'KFold',5,'Stratify',true);
acc = zeros(cvp.NumTestSets,1);

for i=1:cvp.NumTestSets
    tr = training(cvp,i); te = test(cvp,i);

    % Standardize using training set
    mu = mean(X(tr,:),1);
    sigma = std(X(tr,:),[],1);
    Xtr = (X(tr,:) - mu) ./ sigma;
    Xte = (X(te,:) - mu) ./ sigma;
    
    % Train SVM with linear kernel
    Mdl = fitcsvm(Xtr,Y(tr),'KernelFunction','linear');
    Yp  = predict(Mdl,Xte);
    acc(i) = mean(Yp==Y(te));
end
fprintf('5-fold CV accuracy: %.2f%%\n', mean(acc)*100);