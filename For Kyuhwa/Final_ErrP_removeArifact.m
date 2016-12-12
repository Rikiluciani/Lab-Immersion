function [ EEG_clean,EYE_BLINK,EOG,EMG] = Final_ErrP_removeArifact( EEG_raw_errp,fs_errp,chans_errp,F)
%ERRP_ARTIFACTREMOVAL Remove artifact with FORCe assuming that the raw
%data, sampling rate and channel topography is loaded.

windowLength = (0.250*fs_errp); % WINDOW of 250ms
N = windowLength*floor(size(EEG_raw_errp,1)./windowLength); 
                                               % hp: the way to do that is 
                                               % to do an euclidian division 
                                               % between the trial length
                                               % and the window length
                                               % rather than setting it 
                                               % manually
useAcc = 0;
EEG_clean = []; 
EYE_BLINK=zeros(1,floor(size(EEG_raw_errp,1)./windowLength)+1);
EOG=zeros(1,floor(size(EEG_raw_errp,1)./windowLength)+1);
EMG=zeros(1,floor(size(EEG_raw_errp,1)./windowLength)+1);
i = 1;
FLAG=F;
for windowPosition = 1:windowLength:N,
    window = windowPosition:(windowPosition+windowLength)-1;
    % Use FORCe...
    disp([num2str(i), ') ' 'Window from ' num2str(windowPosition/fs_errp),'s to ',num2str((windowPosition+windowLength-1)/fs_errp),'s']);
    tic;
    [cleanEEG,parameters] = Final_FORCe( EEG_raw_errp(window,:)', fs_errp, chans_errp, useAcc );
    disp(['Time taken to clean 250ms EEG = ' num2str(toc) 's.']);
    % Put together the cleaned EEG time series.
    EEG_clean = [EEG_clean cleanEEG];   % Note1: should be same format as input (data x channel)
                                        % why is it transposed ? Related: Note2.
if (F==0)
    FLAG=2;
end
if ((((sum(parameters(6,:))>2 && sum(parameters(7,:))>0)))&&FLAG==2)
        disp('In this 250ms EEG there is EMG artifact'); 
    EMG(1,round(windowPosition/windowLength)+1)=EMG(1,round(windowPosition/windowLength)+1)+1;  
end
if (F==0)
    FLAG=1;
end
if ((sum(parameters(9,:))>0 || sum(parameters(10,:)>0))&&FLAG==1)
        disp('In this 250ms EEG there is eye blink artifact');    
    EYE_BLINK(1,round(windowPosition/windowLength)+1)=EYE_BLINK(1,round(windowPosition/windowLength)+1)+1;
end
if (F==0)
    FLAG=3;
end
if (((sum(parameters(1,:))>0)  && (sum(parameters(11,:))> 0))&&FLAG==3)
        disp('In this 250ms EEG there is EOG artifact'); 
    EOG(1,round(windowPosition/windowLength)+1)=EOG(1,round(windowPosition/windowLength)+1)+1;  
end
i = i+1;

end
EYE_BLINK(end) = EYE_BLINK(end-1);
EOG(end) = EOG(end-1);
EMG(end) = EMG(end-1);
EEG_clean = EEG_clean'; %Note2: transposing so that it has the same shape as the input

end

