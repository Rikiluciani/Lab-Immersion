%% Demo to asses Force with our data

clear all
clc
close all

load('Data_Go.mat');
load('EMG_GT.mat');
load('BLINK_GT.mat');
load('EOG_GT.mat');
SitCalmGT=zeros(15,17);

fs = sample_rate;   clear sample_rate        % Sample rate (512Hz)
fs_errp = fs;

eyeroll_artifact(:,2069,:) = 0;
who('SitCalm_artifact','blink_artifact','eyeroll_artifact','neckLR_artifact');
DATA=input('Select an artifact to process: ');

if(DATA==eyeroll_artifact)
    F=3;
elseif(DATA==neckLR_artifact)
    F=2;
elseif(DATA==blink_artifact)
    F=1;
else
    F=0;
end

B=DATA(64:65,:,:);                        % mean of A1,A2
REF=mean(B);

chan_index=[12 19 20 21 22 23 28 29 30 31 32 37 38 39 40 41];
EEG_raw = DATA(chan_index,:,:);
chans_errp = new_chans(1,chan_index);

t=0:1/fs:(size(EEG_raw,2)-1)/fs;             % t-axis

% Removing of the ear channels A1,A2
for i=1:size(EEG_raw,3)
 for j=1:size(EEG_raw,1)   
    EEG_raw(j,:,i)=EEG_raw(j,:,i)-REF(1,:,i);
 end  
end

%% Plotting the non-filtered EEG data 
%x=input('Select a trial to process (from 1 to 15): ');
for x=1:15
EEG_raw_errp = reshape(EEG_raw(:,:,x),size(EEG_raw,1),size(EEG_raw,2))';

disp( 'Plotting raw data...' ); % ONLY the first 62 channel
figure%(1);
clf;
for chNo = 1:size(EEG_raw_errp,2)

    plot(t,EEG_raw_errp(:,chNo)-((chNo-1)*100));
   
    hold on;
end
title(['EEG and no artifact removal of the trial ' num2str(x)]);
hold off
set(gca,'XGrid','on','Xtick',0:0.25:t(end))

% pause

%% Testing FORCe algorithm

[output(:,:,1),EYE_BLINK,EOG,EMG] = Final_ErrP_removeArifact(EEG_raw_errp,fs,chans_errp,F);

disp('Done. Get back to work.');
asse=0.0:0.25:(floor(size(EEG_raw_errp,1)/(fs*0.25)))*0.25;

if (F==3)
% EOG ARTIFACT
GT = EOG_GT(x,:);
TP = ((GT == 1) & (EOG == 1));
TN = ((GT == 0) & (EOG == 0));
FP = ((GT == 0) & (EOG == 1));
FN = ((GT == 1) & (EOG == 0));
elseif (F==2)
% EMG ARTIFACT 
GT = EMG_GT(x,:);
TP = ((GT == 1) & (EMG == 1));
TN = ((GT == 0) & (EMG == 0));
FP = ((GT == 0) & (EMG == 1));
FN = ((GT == 1) & (EMG == 0));
elseif(F==1)
% EYE BLINK
GT = BLINK_GT(x,:);
TP = ((GT == 1) & (EYE_BLINK == 1));
TN = ((GT == 0) & (EYE_BLINK == 0));
FP = ((GT == 0) & (EYE_BLINK == 1));
FN = ((GT == 1) & (EYE_BLINK == 0));
elseif(F==0)
GT = SitCalmGT(x,:);
TP = ((GT == 1) & (EYE_BLINK == 1)) + ((GT == 1) & (EMG == 1)) + ((GT == 1) & (EOG == 1));
TN = ((GT == 0) & (EYE_BLINK == 0)) + ((GT == 0) & (EMG == 0)) + ((GT == 0) & (EOG == 0));
FP = ((GT == 0) & (EYE_BLINK == 1)) + ((GT == 0) & (EMG == 1)) + ((GT == 0) & (EOG == 1));
FN = ((GT == 1) & (EYE_BLINK == 0)) + ((GT == 1) & (EMG == 0)) + ((GT == 1) & (EOG == 0));
end

hold on
h(1)=stairs(asse,150*GT,'k','LineWidth',2); 
h(2)=stairs(asse,100*TP,'r','LineWidth',2); 
h(3)=stairs(asse,100*FP,'g','LineWidth',2); 
h(4)=stairs(asse,100*FN,'c','LineWidth',2);
h(5)=stairs(asse,100*zeros(1,17),'m','LineWidth',2);

M={'Ground Truth','True Positive','False Positive', 'False Negative'};
legend(h,M,'Location','SouthOutside','Orientation','Horizontal');

% Computing rates
TPR = TP / (TP + FN);
TNR = TN / (FP + TN);
FPR = FP / (FP + TN);
FNR = FN / (TP + FN);
% % % Rates = struct('TPR', TPR, 'TNR', TNR, 'FPR', FPR, 'FNR', FNR); clear TPR TNR FPR FNR
ACC = (TP+TN)/(TP+TN+FP+FN);
title(['EEG and no artifact removal of the trial ' num2str(x) '; ACCURACY: ' num2str(ACC)]);
end
%%
% %% Plotting the filtered data

% % output=output';

% figure(2);
% clf;
% hold on
% t1=0:1/fs:(size(output,1)-1)/fs;             % t1-axis
% for chNo = 1:size(output,2)
%     plot(t1,output(:,chNo)-((chNo-1)*100),'r');
% %    hold on;
% end
% title(['EEG with artifact removed of the trial' num2str(x)]);
% hold off;


