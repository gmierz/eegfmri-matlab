function [EEGres,timep] = remove_gradient(EEG, chanN)
%[EEGres,timep] = remove_gradient(EEG, chanN)
% Removes gradient artifacts from EEG data.
%
%   INPUT:
%       EEG             -  The raw EEG data, should have a 5000Hz sampling rate.
%   OPTIONAL INPUT:
%       chanN           -  A reference channel to use in gradient artifact
%                          subtraction (defaults to channel 48).
%   RETURN:
%       EEGres          -  The cleaned EEG data set.
%       timep           -  Time points of the gradients.
%
% EXAMPLE:
%     [cleaneeg, ~] = remove_gradient(dirtyeeg, 48)
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

    % step 1: get the indices of the gradient epochs 
    hp = eegfiltfft(EEG.data(chanN,:),5000,20,5000);
    allhp = eegfiltfft(EEG.data([1:31],:),5000,30,5000);
    stdhp = std(allhp,0,1);
    dhp = diff(hp,2,2);
    xc = xcorr(hp,5000,'coeff');
    wsize = median(diff(find(xc>.98)));
    allsi = zeros(round(length(hp)./wsize),wsize) ; testepochs = zeros(size(allsi));

    icount = 1;
    for i=1:wsize:length(hp)-wsize*2
        epochi = stdhp(i:i+wsize-1);
        maxi = 100;
        allmaxes(icount) = maxi+i-1;
        allsi(icount,:) = i+maxi:i+maxi+wsize-1;
        testepochs(icount,:) = hp(allsi(icount,:));
        icount = icount + 1;
    end
    zmaxes = zeros(1,size(EEG.data,2)) ; zmaxes(allmaxes) = 1;

    % step 2 find the indices where no gradient is present and remove those:
    sumsv = sum(diff(testepochs,1,2).^2,2);
    medsv = median(sumsv) ; stdsv = std(sumsv); 
    bads = (sumsv < medsv - stdsv | sumsv > medsv + stdsv);
    allsi(bads,:) = []; testepochs(bads,:) = []; teps = testepochs(:,1:10:end);
    corrs = corr(teps'); bad2s = find(zscore(sum(corrs))<-6);
    allsi(bad2s,:) = []; testepochs(bad2s,:) = [];
    timep(1) = allsi(1,end); timep(2) = allsi(end,1);

    % step 3 subtract the average artifact from each channel separately: 
    cleanchans = zeros(size(EEG.data));
    for chanN=1:64 ; disp(['chan=',num2str(chanN)]);
        chan = EEG.data(chanN,:);
        fftchan = fft(chan); freqbinsz = 5000/length(fftchan);
        high_cutoff = 5;
        lowp_ind = round(high_cutoff/freqbinsz);
        lpfft = fftchan;
        lpfft(lowp_ind:end-lowp_ind) = 0;
        lowpower = real(ifft(lpfft));

        hpfft = fftchan;
        hpfft([1:lowp_ind,end-lowp_ind:end]) = 0;
        highpower = real(ifft(hpfft));

        % for all gradient indices, get the artifact at that index:
        allhp = zeros(size(allsi)); allraw = zeros(size(allsi));
        for i=1:size(allsi,1)
            allraw(i,:) = highpower(allsi(i,:));
        end

        filt = fspecial('gaussian',[160,1],80);
        padchan = zeros(size(allraw,1)+1000,size(allraw,2));
        padchan(501:end-500,:) = allraw;
        padchan(1:500,:) = (flipud(allraw(1:500,:)));
        padchan(end-499:end,:) = (flipud(allraw(end-499:end,:)));
        filtchan = imfilter(padchan,filt,'replicate');
        filtchan = filtchan(501:end-500,:);

        subts = zeros(size(chan));
        for i=1:size(allsi,1); subts(allsi(i,:)) = allraw(i,:) - filtchan(i,:); end
        subts_final = subts + lowpower;

        cleanchans(chanN,:) = subts_final;
    end

    EEG2 = EEG;
    EEG2.data = cleanchans; 
    EEGres = pop_resample(EEG2,250);
    timep = timep./20;

end