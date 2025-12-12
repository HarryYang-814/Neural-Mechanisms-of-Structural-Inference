%% ------------------------------------------------------------------------
% SECTION 6 · HELPER FUNCTIONS                                            
% -------------------------------------------------------------------------

function miChan = pac_MI(ep, fs, phBand, amBand, latSmp)
% Compute Tort modulation‑index for each channel in epoch ep
% ep      : [chan × time]
% fs      : sampling‑rate
% phBand  : [low high] Hz for phase
% amBand  : [low high] Hz for amplitude
% latSmp  : samples defining latency window

    nChan = size(ep,1);
    nBins = 36;                           % phase bins (10° width)
    phaseEdges = linspace(-pi, pi, nBins+1);

    % Filter design (2‑pass zero‑phase)
    bpPh = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',phBand(1), ...
                      'HalfPowerFrequency2',phBand(2),'SampleRate',fs);
    bpAm = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',amBand(1), ...
                      'HalfPowerFrequency2',amBand(2),'SampleRate',fs);

    miChan = zeros(1,nChan);
    for ch = 1:nChan
        sig = ep(ch, :);
        ph   = angle(hilbert(filtfilt(bpPh, sig)));
        amp  = abs(  hilbert(filtfilt(bpAm, sig)) );
        ph   = ph(latSmp);
        amp  = amp(latSmp);

        % Bin amplitude by phase
        meanAmp = zeros(1,nBins);
        for b = 1:nBins
            idx = ph >= phaseEdges(b) & ph < phaseEdges(b+1);
            meanAmp(b) = mean(amp(idx));
        end
        meanAmp(isnan(meanAmp)) = 0;
        meanAmp = meanAmp / sum(meanAmp); % probability distribution

        % Modulation Index (KL divergence)
        miChan(ch) = (log(nBins) + sum(fillmissing(meanAmp .* log(meanAmp),'constant',0))) / log(nBins);
    end
end

