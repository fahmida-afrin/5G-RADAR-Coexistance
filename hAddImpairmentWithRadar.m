%% ========================================================================
% hAddImpairmentWithRadar.m - ACTUAL RADAR INJECTION INTO RECEIVED SIGNAL
% 
% THIS IS THE MISSING PIECE!
% 
% This function:
% ✓ Takes transmitted waveform
% ✓ Applies CDL-C channel
% ✓ ADDS RADAR WAVEFORM to received signal
% ✓ Adds AWGN noise
% ✓ Returns degraded signal
%
% This is called for EVERY packet transmission in the network simulator
% ========================================================================

function [rxWaveform, pathGains] = hAddImpairmentWithRadar(txWaveform, varargin)
    % Add channel impairment including radar interference
    %
    % Usage:
    %   [rxWaveform, pathGains] = hAddImpairmentWithRadar(txWaveform)
    %
    % Or if you need access to radar signal:
    %   [rxWaveform, pathGains] = hAddImpairmentWithRadar(txWaveform, ...
    %       radarWaveform, radarFreqOffset, sampleRate)
    
    % Parse inputs
    if nargin >= 4
        radarWaveform = varargin{1};
        radarFreqOffset = varargin{2};
        sampleRate = varargin{3};
    else
        % Default: no radar
        radarWaveform = [];
        radarFreqOffset = 0;
        sampleRate = 30.72e6;  % Standard NR sample rate
    end
    
    % Ensure column vector
    if size(txWaveform, 2) > 1
        txWaveform = txWaveform.';
    end
    
    numSamples = length(txWaveform);
    
    % ====================================================================
    % STEP 1: CDL-C MULTIPATH CHANNEL
    % ====================================================================
    
    % Create impulse response (CDL-C profile)
    % Delays and gains for 6 dominant paths
    pathDelays = [0; 300; 710; 1090; 1730; 2510] * 1e-9;  % seconds
    pathGains_dB = [0; -2.2; -4; -6; -8.2; -9.7];  % dB
    
    % Convert to linear and normalize
    pathGains = 10.^(pathGains_dB / 20);
    pathGains = pathGains / norm(pathGains) * sqrt(0.5);
    
    % Build impulse response
    maxDelay = max(pathDelays);
    numDelayTaps = ceil(maxDelay * sampleRate) + 1;
    h = zeros(numDelayTaps, 1);
    
    for p = 1:length(pathDelays)
        idx = round(pathDelays(p) * sampleRate) + 1;
        if idx <= length(h)
            h(idx) = pathGains(p);
        end
    end
    
    % Apply CDL-C convolution to received signal
    rxWaveform_cdl = conv(txWaveform, h, 'same');
    
    % Normalize to input power
    txPower = mean(abs(txWaveform).^2);
    rxPower = mean(abs(rxWaveform_cdl).^2);
    if rxPower > 0
        rxWaveform_cdl = rxWaveform_cdl * sqrt(txPower / rxPower);
    end
    
    % ====================================================================
    % STEP 2: RADAR WAVEFORM INJECTION (THE CRITICAL STEP!)
    % ====================================================================
    
    rxWaveform_interference = rxWaveform_cdl;
    
    if ~isempty(radarWaveform) && length(radarWaveform) > 0
        % Ensure radar is column vector
        if size(radarWaveform, 2) > 1
            radarWaveform = radarWaveform.';
        end
        
        % Truncate or repeat radar to match tx length
        if length(radarWaveform) >= numSamples
            radarSeg = radarWaveform(1:numSamples);
        else
            % Repeat radar signal to fill length
            nReps = ceil(numSamples / length(radarWaveform));
            radarRep = repmat(radarWaveform, nReps, 1);
            radarSeg = radarRep(1:numSamples);
        end
        
        % Apply frequency shift to radar signal
        % Shift by radarFreqOffset (typically 3 MHz)
        t = (0:numSamples-1).' / sampleRate;
        freqShift = exp(1j * 2 * pi * radarFreqOffset * t);
        radarShifted = radarSeg .* freqShift;
        
        % Scale radar signal to received power
        % Radar power after path loss (assume 5 km distance)
        radarTxPower = 1000;  % 1 kW
        distance = 5000;  % 5 km
        c = 3e8;
        freq = 3.595e9 + radarFreqOffset;  % Radar center frequency
        wavelength = c / freq;
        pathLoss = (4 * pi * distance / wavelength)^2;
        radarRxPower = radarTxPower / pathLoss;
        
        % Scale to match radar received power
        radarPower = mean(abs(radarShifted).^2);
        if radarPower > 0
            radarScaled = radarShifted * sqrt(radarRxPower / radarPower);
        else
            radarScaled = radarShifted;
        end
        
        % THIS IS THE CRITICAL LINE:
        % ADD RADAR DIRECTLY TO RECEIVED SIGNAL
        rxWaveform_interference = rxWaveform_cdl + radarScaled;
        
        fprintf('DEBUG: Radar injected. Radar power: %.2e, Signal power: %.2e\n', ...
            radarRxPower, mean(abs(rxWaveform_cdl).^2));
    end
    
    % ====================================================================
    % STEP 3: AWGN NOISE
    % ====================================================================
    
    % SNR = 15 dB (typical operating point)
    SNRdB = 15;
    
    % Calculate noise power
    signalPower = mean(abs(rxWaveform_interference).^2);
    SNRLinear = 10^(SNRdB/10);
    noisePower = signalPower / SNRLinear;
    
    % Generate complex Gaussian noise
    noise = (randn(numSamples, 1) + 1j*randn(numSamples, 1)) / sqrt(2);
    noise = noise * sqrt(noisePower);
    
    % ADD NOISE TO SIGNAL
    rxWaveform = rxWaveform_interference + noise;
    
end

%% ========================================================================
% THIS FUNCTION IS CALLED BY NETWORK SIMULATOR FOR EVERY TRANSMISSION
% ========================================================================
%
% The network simulator automatically calls this for:
% - Every DL packet (gNB → UE)
% - Every UL packet (UE → gNB)
%
% Each call:
% 1. Receives txWaveform (what was transmitted)
% 2. Applies CDL-C fading
% 3. ADDS RADAR to the signal  <-- THIS IS THE KEY LINE
% 4. Adds noise
% 5. Returns rxWaveform (what was received)
%
% The received waveform is then passed to:
% - OFDM demodulator
% - Equalization
% - CQI calculation
% - BLER/PER determination
% - HARQ decision
%
% So the radar degradation AUTOMATICALLY affects:
% ✓ SINR (measurement from received signal)
% ✓ Decoding (more errors due to interference)
% ✓ HARQ (more retransmissions)
%
% For BOTH UL and DL equally!

