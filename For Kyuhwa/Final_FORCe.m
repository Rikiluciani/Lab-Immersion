function [cleanEEG,parameters] = Final_FORCe( EEGdata, Fs, chanLocs, useAcc )
    %
    % FORCe
    %
    %   Fully automated and Online artifact Removal for brain-Computer
    %   interfacing.
    %
    % First perform hard threhsolding to reject very high amplitude
    % channels. Then re-reference, calculate wavelet coefficients,
    % estimate (via SOBI) ICA de-mixing matrix, extract features, apply
    % thresholds. Flag ICs for inclusion or exclusion, soft threshold
    % spikes and reconstruct cleaned signals.
    %
    % Inputs:
    %
    %   EEGdata     - M x N matrix of M channels by N samples. Note, the
    %           method is currently hard coded to deal with 1 s of EEG.
    %   accSigs     - Signals from the accelerometer (currently not fully
    %           tests; the acceleromter didn't work when attempting to
    %           acquire test signals).
    %   Fs          - sampling rate (note, only tested with 512 Hz so far,
    %           some parts may expect this).
    %   chanLocs    - Locations of the channels used in EEGlab format.
    %   useAcc      - 0 = don't use, 1 = use (I suggest setting to zero
    %           until someone has implemented and tested this part).
    %
    % Output:
    %
    %   cleanEEG    - M x N matrix of cleaned EEG where N is of length 100
    %      ms, only the most recent 100 ms are cleaned.
    %
    % Author: Ian Daly, 2013
    %
    % License:
    %
    % This program is free software; you can redistribute it and/or modify it under
    % the terms of the GNU General Public License as published by the Free Software
    % Foundation; either version 2 of the License, or (at your option) any later
    % version.
    %
    % This program is distributed in the hope that it will be useful, but WITHOUT
    % ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    % FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License along with
    % this program; if not, write to the Free Software Foundation, Inc., 59 Temple
    % Place - Suite 330, Boston, MA  02111-1307, USA.
    %
    %
    % REFERENCES:
    %
    % Please cite the following manuscript when using this method.
    %
    % Daly, I. et al., 2014. FORCe: Fully Online and automated artifact 
    %   Removal for brain-Computer interfacing. IEEE transactions on neural
    %   systems and rehabilitation engineering : a publication of the IEEE 
    %   Engineering in Medicine and Biology Society. Available at: 
    %   http://www.ncbi.nlm.nih.gov/pubmed/25134085
    %
    %
    %**********************************************************************
       
    % Begin function proper, check the dimmensions of the input data.   
    N = size( EEGdata,1 ); %number of channels
    M = size( EEGdata,2 ); %number of samples
    
    % Check for NaNs and Infs, if found display warning.
    for chNo = 1:N,
        if any( isnan( EEGdata(chNo,:) ) ),
            disp( 'Warning: signal contains NaNs' );
        end
        if any( isinf( EEGdata(chNo,:) ) ),
            disp( 'Warning: signal contains Infs' );
        end
    end
    
    % Step 1. Hard threshold channels for extreme values (> 200uV).
    remCh = [];
    for chNo = 1:N,
        if max( EEGdata(chNo,:) ) > 200,
            remCh = [remCh chNo];
        end
    end
    
    
    % Check we have some usable data (ie. there are some channels with
    % amplitude below the 200uV threshold). If no usable date found just
    % retunr what we have (ie. retunr input without cleaning).
    if length( unique( remCh ) ) == size( EEGdata,1 ),
        disp( 'Error: No usable data in this EEG epoch: aborting!' );
        cleanEEG = EEGdata;
        return;
    end
    
    % Identify channels in frontal locations ...
    frontChans = [];
    for chNo = 1:length( chanLocs ),
        if length( strfind( chanLocs(chNo).labels, 'F' ) ) > 0,
            frontChans = [frontChans chNo];
        end
    end

  
    % ... and other locations.
 
    % otherChans = setdiff( 1:16,frontChans ); %%
    NCh = length( chanLocs );
    otherChans = setdiff( 1:NCh,frontChans );
   
    % Estimate the new scalp data from the old minus removed channels.
    chData = estimateRemovedChs( EEGdata,remCh,chanLocs );
    
    
    % START OF WAVELET DECOMPOSITION AND INDEPENDENT COMPONENT ANALYSIS
    % Estimate wavelet coefficients and ICs via Smlet wavelets and SOBI.
    [ICs mixMat wavePacks] = applyOnlineICAmethodWaveNowsobi( chData );
    wavePacksOut = wavePacks;
    tNodes = get( wavePacks{1},'tN' );
    
    tUse = 1:length( tNodes );
    
    
    % ALL STARTS FROM HERE
    % Step through the terminal nodes.
    for tN = 1:length( tUse ),
        
        % Step 4(a). Extract features relevent to acceleromter measures.
        if useAcc == 1,
            [featsAcc sigsAcc] = getFeaturesAcc( ICs{tUse(tN)},accSigs,Fs );
        else
            featsAcc = [];
            sigsAcc = [];
        end
        
        if tN == 1,% If this is the first terminal node (the approximation coefficients).
            
            % Extract features to be thresholded.
            clear temp;
            remICsPF = []; % if the amplitude of the projection exceed 100uV.
            remICsP2P = [];% if the peak to peak distance between maxima and
                            % minima amplitudes is major than 60.
            remICsEOG = [];
            
            % For each IC (estimated from SOBI ICA).
            for iNo = 1:size( ICs{tUse(tN)},1 ),
                % Estimate the IC projection onto the scalp.
                ICtemp = ICs{tUse(tN)};
                for iT = 1:size( ICtemp,1 ),
                    if iT ~= iNo,
                        ICtemp(iT,:) = zeros(1,size(ICtemp,2));
                    end
                end
                projIC{iNo} = mixMat{tUse(tN)} * ICtemp;
                
                % Check the projections of each IC against some
                % simple thresholds.
                for pNo = 1:size( projIC{iNo},1 ),
                    % check thresholds, does the amplitude of the
                    % projection exceed 100uV.
                    passFail(pNo) = checkThresholds(projIC{iNo}(pNo,:));
                    % Check the peak to peak distance between maxima and
                    % minima amplitudes.
                    p2p(pNo) = max( projIC{iNo}(pNo,:) ) - min( projIC{iNo}(pNo,:) );
                end
                % If any amplitudes exceed this threshold mark the
                % corresponding IC fro removal.
                if any( passFail == 0 ),
                    remICsPF = [remICsPF iNo];
                end
                % Same for p2p
                if any( p2p > 90 ), %try to change this value
                    remICsP2P = [remICsP2P iNo];
                end
                if any( p2p > 56 & p2p < 90)
                %if any( p2p > 50 & p2p < 80),% && any (p2p < 50 ), %try to change this value
                    remICsEOG = [remICsEOG iNo];
                end
               
                % Get kurtosis values over projections
                ent(iNo) = mean( kurtosis( projIC{iNo}' ) ); %kurtosis of the scalp projection
                                                            %of a single IC
                                                            
                % Check PSDs.
                clear ps;
                clear ps2;
                for pNo = 1:size( projIC{iNo},1 ),
                    % Get power spectrum OF THE REAL IC
                    psT = powerspectrum( projIC{iNo}(pNo,:)',Fs );
                    
                    % Match to 1/f distribution. IDEAL 1/F DISTRIBUTION
                    clear idealDistro;
                    idealDistro(1,:) = psT(1,:); %the same frequency axis
                    idealDistro(2,2:end) = 1./psT(1,2:end);
                    % Match to magnitude of other distros (NORMALISE BOTH).
                    idealDistro(2,2:end) = idealDistro(2,2:end) ./ max( idealDistro(2,2:end) );
                    psCheck = psT(2,2:end) ./ max( psT(2,2:end) );
                    %figure, plot(idealDistro(1,2:end),psCheck,'r'), hold on,plot(idealDistro(1,:),idealDistro(2,:),'g');
                    
                    % Calc. euclidean distance between ideal (1/F) spectra
                    % and measured spectra.
                    specDistT(pNo) = dist( psCheck,idealDistro(2,2:end)' );
                    
                    %Calculate the mean PSD of the scalp projection of each
                    %IC in frequencies above 30Hz
                    gammaPSDT(pNo) = mean( psT(2,find(psT(1,:) > 30 )) );
                    slow_wave(pNo) = mean( psT(2,find(psT(1,:) < 9  )) );
                end
                % Take mean distance over projects.
                specDist(iNo) = mean(specDistT);
                % Maximum PSD in gamma and above bands.
                [gammaPSD(iNo),~] = max( gammaPSDT );
                % Maximum PSD in lower bands
                [lowerPSD(iNo),~] = max( slow_wave);
                
                % Check stds of projections of ICs.
                for pNo = 1:size( projIC{iNo},1 ),
                    stdProjT(pNo) = std( projIC{iNo}(pNo,:) );
                end
                stdProj(iNo) = max( stdProjT );
                
                % Check std and std ratio (frontal channels / other channels)
                for pNo = 1:size( projIC{iNo},1 ),
                    stdValueT(pNo) = std( projIC{iNo}(pNo,:) );
                end
                stdRatio(iNo) = mean( stdValueT(frontChans) ) / mean( stdValueT(otherChans ) );
               
                % Calculate AMI of each IC. AMI=Auto Mutual Information
                featsClust(iNo,:) = extractFeaturesMultiChsWaveAMI( projIC{iNo},Fs );
                
            end
               
             % Threshold Std ratio. Larger standard deviation of amplitude 
             % on frontal channels may indicate the presence of EOG activity
             % -> so we calculate R(stdRatio)
             stdRatioRem = find( stdRatio > (mean(stdRatio)+(0.65*std(stdRatio))) ); %0.65 instead of 1 
                          
             % Threshold Gamma rem. High power PSD in the gamma frequency
             % band and above could indicate the presence of EMG artifact
             % contamination in the EEG.
             gammaRem = find( gammaPSD > 2.4); %1.7 );
             
             lowerRem = find( lowerPSD > 5.0);
             
             % Check spikeness of data.  LOOK WELL AT THIS PART
             if tN == 1,
                 for iNo = 1:size( ICs{tN},1 ),
                     % Find spike zones in data.
                     muSig = abs( mean( projIC{tN} ) );
                     A1 = find( muSig(2:end-1) > muSig(3:end) )+1;
                     spikePos = A1(find( muSig(A1) > muSig(A1-1) ));
                     % Calculate coefficients of variation for each spike zone.
                     coefVar = zeros(1,length(spikePos));
                     for i = 1:length( spikePos ),
                         coefVar(i) = std(abs(mean(projIC{tN}(:,spikePos(i)-1:spikePos(i)+1)))) / mean(abs(mean(projIC{tN}      (iNo,spikePos(i)-1:spikePos(i)+1))));
                     end
                     % Apply soft thresholding to any spike zone for which the
                     % coefficient of variation exceeds the threshold.
                     T = 0.1 * (mean(coefVar)+std(coefVar));
                     if length( coefVar ) > 0,
                         noSpikes(iNo) = length( find( coefVar > T ) ) / length( coefVar );
                     else
                         noSpikes(iNo) = 0;
                     end
                 end
                 remNoSpikes = find( noSpikes >= 0.25 ); %it is the second way to find spiking activities
              end
             
             % Threshold the kurtosis values  
             remICsKurt = [];  % here the IC's projection that not respect the condition.
             for i = 1:length( ent ),
                if ent(i) > mean((ent)) + (0.5*std((ent))) || ent(i) < mean(ent) - (0.5*std(ent)),
                    remICsKurt = [remICsKurt i];
                end
             end
             
             % Threshold the distance to the 1/F distribution. 
             rems1F = find( specDist > 3.5); 
             %the power spectra of the clean EEG is known to approximately
             %follow a 1/F distribution
             %A very wide threshold is applied (it is an approximate rule).
                        
             
             % Theshold the ratio of PSDs between < 20Hz and > 20Hz.????
             for i = 1:size( ICs{tN},1 ),
                ps = powerspectrum( ICs{tN}(i,:)',Fs );
                muLow = mean(ps(2,find( ps(1,:) < 20 )));
                muHigh = mean(ps(2,find( ps(1,:) > 20 )));
                psRatio(i) = muHigh / muLow;
            end
            remICspsRatio = find( psRatio > 1.0 );
                        
            % Threshold the std of the IC projects.
            remSTD = find( stdProj > (mean(stdProj)+(1.65*std(stdProj))) );
            %High standard deviation in the EEG have also been reported to
            %indicate the presence of EMG and other artifacts.
            
            
            % Threshold the spikiness (std) of the ICs.
            % Spiking activity in the signal may indicate the presence of
            % EMG artifacts in an IC -> 2 ways to identify them.
            remSpike = [];
            for i = 1:size( ICs{tN},1 ),
                if max( abs( ICs{tN}(i,:) ) ) > (mean( abs(ICs{tN}(i,:)) ) + (3*std( (ICs{tN}(i,:)) )))
                    remSpike = [remSpike i];  % first way for spiking activity
                end
            end
            
            % Threshold the AMI between 2.0 ... 3.0 
            %They are the 2 thresholds  for the AMI parameter
            remICCs = [find( featsClust > 3 ); find(featsClust < 2.0)]';
             
            % Check which ICs to remove. Remove ICs for which more than 3
            % thresholds are exceeded.
            remICs = ( [gammaRem remICsKurt remICCs remSTD stdRatioRem remICsPF remICsP2P remNoSpikes remICspsRatio remSpike rems1F remICsEOG lowerRem] ); %  remPSDr  remICCs remSpike  rems1F
           
            % Vote.
            
            % Now I search how many ICs to remove there are.
            %for iNoType = 1:16,
            for iNoType = 1:NCh,
                lenINo(iNoType) = length( find( remICs == iNoType ) );
            end
            remICs = find( lenINo > 3 );
            
            parameters=zeros(length(remICs),11);
            j=1;
            while j<=length(remICs)
                i=remICs(j);
                %parameters(j,1)=sum(remICCs==i)>0; % AMI 
                parameters(j,1)=sum(lowerRem==i)>0;
                
                parameters(j,2)=sum(remNoSpikes==i)>0; %NoSpikes 2a spiking activities
                     
                parameters(j,3)= sum(remSpike==i)>0; %Spike 2b spiking activities
                
                parameters(j,4)= sum(remICsKurt==i)>0; %Kurtosis
                
                parameters(j,5)=sum(rems1F==i)>0;  %1/F PSD
                
                parameters(j,6)= sum(gammaRem==i)>0; % Gamma PSD
                
                parameters(j,7)= sum(remSTD==i)>0;  % STD
                   
                parameters(j,8)=sum(stdRatioRem==i)>0;% STD RATIO
                    
                parameters(j,9)=sum(remICsPF==i)>0; %PF peak amplitudes
                    
                parameters(j,10)=sum(remICsP2P==i)>0; %P2p peak amplitudes
                                    
                %parameters(j,11)=sum(remICspsRatio==i)>0; % ??? 
                parameters(j,11)=sum(remICsEOG==i)>0;
                j=j+1;
            end
            parameters=parameters';
            
            % Check number of removed ICs is not the same as the total number
            % of ICs.
            if length( remICs ) == size( ICs{tUse(tN)},1 ),
                disp( 'Warning: All ICs removed!' );
            end
            
            % Remove the ICs
            keepICs = ICs{tUse(tN)};
            keepICs(remICs,:) = zeros(length(remICs),size(ICs{tUse(tN)},2));
            a = mixMat{tUse(tN)};
            
        end
       
        % SPIKE ZONE THRESHOLDING -> Within both the approximation and
        % detail coefficients derived from the wavelet decomposition, a
        % soft thresholding approach is applied to attempt to reduce the
        % influence of EMG artifact contamination in the EEG.
        
        % If we're currently dealing with the approximation coefficients
        % then remove the ICs flagged by the thresholding above.
        if tN == 1,
            newSigNode{tUse(tN)} = a * keepICs;
        else
            newSigNode{tUse(tN)} = ICs{tUse(tN)};
        end
        
        % Apply soft thresholding to adjust spike magnitudes.
        if tN == 2 || tN == 3 || tN == 4,
            % Thresholds for detail coefficients.
            checkVal = 0.2;
            adjustVal = 0.07;
        else
            % Thresholds for approximation coefficients.
            checkVal = 0.7;
            adjustVal = 0.8;
        end
        allCoefVar = [];
        clear spikePos;
        clear coefVar;
        for iNo = 1:size( newSigNode{tUse(tN)},1 ),
            
            muSig = newSigNode{tUse(tN)}(iNo,:);
            A1up = find( muSig(2:end-1) > muSig(3:end) )+1;
            A1lo = find( muSig(2:end-1) < muSig(3:end) )+1;
            spikePos{iNo} = [A1up(find( muSig(A1up) > muSig(A1up-1) )) ...
                A1lo(find( muSig(A1lo) < muSig(A1lo-1) ))];
            indUp = find(spikePos{iNo} > 3 );
            indUse = indUp(find( spikePos{iNo}(indUp) < (length(muSig)-2)));
            spikePos{iNo} = spikePos{iNo}( indUse );
            
            % Use wide segments for detection of spike zone coefficients
            % (specifically used +- 3 samples).
            upperVals = var( [newSigNode{tUse(tN)}(iNo,spikePos{iNo}-3); ...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}-2);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}-1);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo});...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}+1);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}+2);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}+3)] );
            lowerVals = std( [newSigNode{tUse(tN)}(iNo,spikePos{iNo}-3); ...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}-2);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}-1);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo});...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}+1);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}+2);...
                newSigNode{tUse(tN)}(iNo,spikePos{iNo}+3)] );
            
            % Coefficient of variation.
            coefVar{iNo} = upperVals ./ lowerVals;
            
            % Collect...
            allCoefVar = [allCoefVar coefVar{iNo}];
        end
        % For each channel apply soft thresholding.
        for iNo = 1:size( newSigNode{tUse(tN)},1 ),
            T = checkVal * (mean(allCoefVar)+std(allCoefVar));
            for i = 1:length( spikePos{iNo} ),
                if abs(coefVar{iNo}(i)) > T,
                    newSigNode{tUse(tN)}(iNo,spikePos{iNo}(i)) = (adjustVal)*newSigNode{tUse(tN)}(iNo,spikePos{iNo}(i));
                end
            end
        end        
    
    end
    
    % Reconstruct from wavelet decompositions.
    for chNo = 1:size( chData,1 ),
        for tN = 1:length( tUse ),
            wavePacksOut{chNo} = write( wavePacksOut{chNo},'data',tNodes(tUse(tN)),newSigNode{tUse(tN)}(chNo,:) );
        end
    end
    for chNo = 1:size( chData,1 ),
        cleanEEG(chNo,:) = wprec( wavePacksOut{chNo} );
    end

    % The below line can be used to plot the original epoch and the cleaned
    % epoch to check the performance of the method (just add a breakpoint
    % and run it :).
    % clf;t=0:1/Fs:1.0; subplot(2,1,1);hold on;for chNo=1:16,plot(t,EEGdata(chNo,:)-(chNo*50));end;hold off;axis([0 1 -900 0]);subplot(2,1,2);hold on;for chNo=1:16,plot(t,cleanEEG(chNo,:)-(chNo*50));end;hold off;axis([0 1 -900 0]);
    
end
