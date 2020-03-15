function [neweeg] = eegfmri_clean(EEG, output_name, minpeaks, epochls, gradientchan)
%[neweeg] = eegfmri_clean(EEG, output_name, minpeaks, epochls, gradientchan)
% Cleans EEG data that has fMRI gradient and BCG aritfacts.
%
% This function can be used to clean EEG data that is contaminated by
% gradient, and balistocardiogram (BCG) artifacts. It performs gradient
% artifact cleaning first, then it cleans out the BCG artifact. See the
% README.md file for more information about how you should use this
% function. It's mostly automatic, but it requires user input when deciding
% which BCG component we should use to determine where the BCG artifacts
% occur. You should pick one of the subplots that show the "nicest" peaks,
% they can be either maximum's or minimum's, you will set the polarity to
% -1 if the peaks are minimums.
%
%   INPUT:
%       EEG             -  The raw EEG data, should have a 5000Hz sampling rate.
%   OPTIONAL INPUT:
%       output_name     -  The name of the output of the cleaned EEG, if it's
%                          not given, the data won't be saved.
%       minpeaks        -  The minimum distance (in points) between the BCG
%                          artifacts (defaults to 180 points).
%       epochls         -  Size of the BCG artifact to extract (defaults to
%                          [100, 185]). Format is [points_before,
%                          points_after].
%       gradientchan    -  A reference channel to use in gradient artifact
%                          subtraction (defaults to channel 48).
%   RETURN:
%       neweeg          -  The cleaned EEG data set.
%
% EXAMPLE:
%     cleaneeg = eegfmri_clean(dirtyeeg, 'clean_eeg.set')
%
% Author: Russell Butler, & Gregory Mierzwinski, Sherbrooke, QC, 03/19/2019

% Copyright (C) Russell Butler, & Gregory Mierzwinski, gmierz1@live.ca
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    % Setup optional arguments
    if ~exist('output_name', 'var') || isempty(output_name)
        output_name = '';
        disp('Output name was not provided, not saving data');
    end
    if ~exist('minpeaks', 'var') || isempty(minpeaks)
        minpeaks = 180;
    end
    if ~exist('epochls', 'var') || isempty(epochls)
        epochls = [100,200-15];
    end
    if ~exist('gradientchan') || isempty(gradientchan)
        gradientchan = 48;
    end
    
    % Remove the gradient
    [correeg, ~] = remove_gradient(EEG, gradientchan);

    filt = eegfiltfft(correeg.data,correeg.srate,1,128); % high pass
    [weights,sphere] = runica(filt(:,1:4:end),'maxsteps',128) ; % ICA
    allweights = weights*sphere; % save weights
    acts = weights*sphere*filt;

    % visualize top 12 components
    figure,for i=1:12 ; subplot(3,4,i); plot(squeeze(acts(i,10000:40000))); title(i); end;
    EEG = correeg;
    bdat = eegfiltfft(EEG.data,EEG.srate,1,128); % for average subtraction (we are not removing values under 1Hz)
    hdat = eegfiltfft(EEG.data,EEG.srate,3,128); % for individual BCG epoch correlation (more accurate to exclude low-freq)
    
    % BCG ICA component (set after visualizing results in line 24)
    bcg_comp = input('Select your BCG component:');
    % polarity of component (peaks positive (1) or negative (-1)) (set after visualizing)
    polarity = input('Select direction of BCG component (-1 if peaks are minimums, 1 otherwise):');

    acts = squeeze(allweights(:,:))*bdat;
    template = acts(bcg_comp,:)*polarity ; % get the bcg component time series selected in previous script
    [pks,locs] = findpeaks(smooth(template),'MINPEAKDISTANCE',minpeaks); % get the peaks, use subject specific value for minpeaks

    elength = epochls(1) + epochls(2); % length of one bcg epoch (samples)
    eprev = epochls(1); epost = epochls(2); % before and after peak number of samples
    eps = zeros(64,length(locs),elength+1); % lower highpass epochs
    heps = zeros(64,length(locs),elength+1); % higher highpass epochs
    epinds = zeros(length(locs),elength+1); % epoch indices
    for i=2:length(locs)-2 % get the BCG epochs using locs from findpeaks
        eps(:,i,:) = bdat(:,locs(i)-eprev:locs(i)+epost); 
        heps(:,i,:) = hdat(:,locs(i)-eprev:locs(i)+epost); 
        epinds(i,:) = locs(i)-eprev:locs(i)+epost;
    end
    meaneps = squeeze(mean(eps,2)); % mean lower highpass epochs (all channels, mean across epoch)


    [bsv,bsi] = sort(sum(abs(meaneps),2),'descend'); % sort channels by epoch magnitude 
    clear corrs
    for i=1:size(heps,1) % create correlation matrix by correlating all epochs
        corrs(i,:,:) = corr(squeeze(heps(i,:,:))');  

    end

    mcorrs = squeeze(mean(corrs(:,:,:),1)); % mean across epochs, can use bsi to only average channels with more BCG
    [sv,si] = sort(mcorrs,2,'descend') ; % sort according to highest correlation for each epoch
    neweeg = EEG; % create new EEG for subtracted data
    for i=2:size(epinds,1)- 2
        for j=1:64
            neweeg.data(j,epinds(i,:)) = squeeze(EEG.data(j,epinds(i,:))) - ...
                squeeze(mean(eps(j,si(i,5:40),:),2))' ; % subtract mean of top correlating (5-40) epochs
        end
    end

    if ~isempty(output_name)
        pop_saveset(neweeg, 'filename', output_name)
    end
end


