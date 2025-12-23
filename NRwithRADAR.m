%% ========================================================================
% MAIN_ACTUAL_RADAR_INJECTION.m
%

clear; clc; close all;
wirelessnetworkSupportPackageCheck

%% ========================================================================
% SECTION 1: CONFIGURATION
% ========================================================================

fprintf('CONFIGURING SIMULATION...\n\n');

rng('default');
simParameters = [];
simParameters.NumFramesSim = 10;
simParameters.SchedulingType = 1;
simParameters.NumUEs = 4;

% UE POSITIONS
simParameters.UEPosition = [100 0 0; 150 0 0; 300 0 0; 400 0 0];

% 5G PARAMETERS
simParameters.NumRBs = 25;
simParameters.SCS = 15;
simParameters.DLBandwidth = 5e6;
simParameters.ULBandwidth = 5e6;
simParameters.DLCarrierFreq = 3.595e9;
simParameters.ULCarrierFreq = 3.595e9;

fprintf('5G: %.1f MHz, %.1f MHz BW, %d RBs\n', ...
    simParameters.DLCarrierFreq/1e6, simParameters.DLBandwidth/1e6, simParameters.NumRBs);

% TDD PATTERN
simParameters.DLULPeriodicity = 5;
simParameters.NumDLSlots = 2;
simParameters.NumDLSyms = 8;
simParameters.NumULSyms = 4;
simParameters.NumULSlots = 2;

% SCHEDULER
simParameters.SchedulerStrategy = 'PF';
simParameters.TTIGranularity = 4;
simParameters.RBAllocationLimitUL = 15;
simParameters.RBAllocationLimitDL = 15;

% CHANNEL
simParameters.ChannelUpdatePeriodicity = 0.2;
simParameters.CQIDelta = 1;
simParameters.CQIvsDistance = [200 15; 300 12; 500 10; 1000 8; 1200 7];

% DMRS
simParameters.DMRSTypeAPosition = 2;
simParameters.PUSCHDMRSAdditionalPosTypeB = 0;
simParameters.PUSCHDMRSAdditionalPosTypeA = 0;
simParameters.PUSCHDMRSConfigurationType = 1;
simParameters.PDSCHDMRSAdditionalPosTypeB = 0;
simParameters.PDSCHDMRSAdditionalPosTypeA = 0;
simParameters.PDSCHDMRSConfigurationType = 1;

% TRAFFIC
simParameters.dlAppDataRate = [160e3, 30e3, 100e3, 10e3];
simParameters.ulAppDataRate = [160e3, 30e3, 100e3, 10e3];

% VISUALIZATION
simParameters.CQIVisualization = false;
simParameters.RBVisualization = false;
enableTraces = true;
simParameters.NumMetricsSteps = 20;

%% ========================================================================
% SECTION 2: RADAR CONFIGURATION & WAVEFORM GENERATION
% ========================================================================

fprintf('CONFIGURING AND GENERATING RADAR WAVEFORM...\n\n');

simParameters.radar = [];
simParameters.radar.prf = 1e3;
simParameters.radar.FreqOffset = 3e6;
simParameters.radar.PulseWidth = 70e-6;
simParameters.radar.BW = 5e6;
simParameters.radar.Power = 1000;
simParameters.radar.SampleRate = 30.72e6;  % Match NR sample rate

fprintf('Radar: %.1f MHz offset, %.1f MHz BW, %.0f kW power\n\n', ...
    simParameters.radar.FreqOffset/1e6, simParameters.radar.BW/1e6, ...
    simParameters.radar.Power/1000);

% Generate LFM radar waveform
fs_radar = simParameters.radar.SampleRate;
pulseWidth = simParameters.radar.PulseWidth;
bw = simParameters.radar.BW;

numPulseSamples = round(pulseWidth * fs_radar);
t_pulse = (0:numPulseSamples-1) / fs_radar;

% LFM: quadratic phase for frequency sweep
chirpRate = bw / pulseWidth;
phase_lfm = pi * chirpRate * t_pulse.^2;
radar_pulse = exp(1j * phase_lfm);

% Window to reduce sidelobes
window = hamming(numPulseSamples)';
radar_pulse = radar_pulse .* window;
radar_pulse = radar_pulse / max(abs(radar_pulse));

% Create radar pulse train over simulation duration
frameDuration = simParameters.NumFramesSim * 10e-3;  % 100 ms
numFrameSamples = round(frameDuration * fs_radar);
pri = 1 / simParameters.radar.prf;
numSamplesPerPRI = round(pri * fs_radar);

radarWaveformFull = zeros(1, numFrameSamples);
for pulseNum = 0:floor(numFrameSamples / numSamplesPerPRI)
    startIdx = pulseNum * numSamplesPerPRI + 1;
    endIdx = min(startIdx + numPulseSamples - 1, numFrameSamples);
    if startIdx <= numFrameSamples
        len = min(numPulseSamples, endIdx - startIdx + 1);
        radarWaveformFull(startIdx:startIdx+len-1) = radar_pulse(1:len);
    end
end

simParameters.radar.waveform = radarWaveformFull;

fprintf('Radar LFM waveform generated: %d samples\n', numFrameSamples);
fprintf('Ready to inject into signal path\n\n');

%% ========================================================================
% SECTION 3: RLC & LOGICAL CHANNELS
% ========================================================================

fprintf('SETTING UP RLC & LOGICAL CHANNELS...\n\n');

numLogicalChannels = 1;
simParameters.LCHConfig.LCID = 4;

lchInfo = repmat(struct('RNTI',[],'LCID',[],'EntityDir',[]), [simParameters.NumUEs 1]);
for idx = 1:simParameters.NumUEs
    lchInfo(idx).RNTI = idx;
    lchInfo(idx).LCID = simParameters.LCHConfig.LCID;
    lchInfo(idx).EntityDir = 2;
end

simParameters.RLCConfig.EntityType = 2;
rlcChannelConfigStruct.LCGID = 1;
rlcChannelConfigStruct.Priority = 1;
rlcChannelConfigStruct.PBR = 8;
rlcChannelConfigStruct.BSD = 10;
rlcChannelConfigStruct.EntityType = simParameters.RLCConfig.EntityType;
rlcChannelConfigStruct.LogicalChannelID = simParameters.LCHConfig.LCID;
simParameters.maxRLCSDULength = 9000;

% Initialize CQI
maxUECQIs = zeros(simParameters.NumUEs, 1);
for ueIdx = 1:simParameters.NumUEs
    matchingRowIdx = find(simParameters.CQIvsDistance(:, 1) > simParameters.UEPosition(ueIdx, 1));
    if isempty(matchingRowIdx)
        maxUECQIs(ueIdx) = simParameters.CQIvsDistance(end, 2);
    else
        maxUECQIs(ueIdx) = simParameters.CQIvsDistance(matchingRowIdx(1), 2);
    end
end

simParameters.InitialChannelQualityUL = zeros(simParameters.NumUEs, simParameters.NumRBs);
simParameters.InitialChannelQualityDL = zeros(simParameters.NumUEs, simParameters.NumRBs);

for ueIdx = 1:simParameters.NumUEs
    simParameters.InitialChannelQualityUL(ueIdx, :) = randi([1 maxUECQIs(ueIdx)], 1, simParameters.NumRBs);
    simParameters.InitialChannelQualityDL(ueIdx, :) = simParameters.InitialChannelQualityUL(ueIdx, :);
end

fprintf('RLC configured\n');
fprintf('CQI initialized\n\n');

%% ========================================================================
% SECTION 4: CREATE GNB & UEs
% ========================================================================

fprintf('CREATING GNB AND UEs...\n\n');

simParameters.DuplexMode = 1;
simParameters.NCellID = 1;
simParameters.Position = [0 0 0];

% Create gNB
gNB = hNRGNB(simParameters);
scheduler = hNRSchedulerRoundRobin(simParameters);
addScheduler(gNB, scheduler);
gNB.PhyEntity = hNRGNBPhy(simParameters);
configurePhy(gNB, simParameters);
setPhyInterface(gNB);

% Create UEs
UEs = cell(simParameters.NumUEs, 1);
for ueIdx = 1:simParameters.NumUEs
    simParameters.Position = simParameters.UEPosition(ueIdx, :);
    UEs{ueIdx} = hNRUE(simParameters, ueIdx);
    UEs{ueIdx}.PhyEntity = hNRUEPhy(simParameters, ueIdx);
    configurePhy(UEs{ueIdx}, simParameters);
    setPhyInterface(UEs{ueIdx});
    
    % Initialize channel quality
    channelQualityInfoUL = struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityUL(ueIdx, :));
    updateChannelQualityUL(gNB.MACEntity.Scheduler, channelQualityInfoUL);
    
    channelQualityInfoDL = struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityDL(ueIdx, :));
    updateChannelQualityDL(gNB.MACEntity.Scheduler, channelQualityInfoDL);
    updateChannelQualityDL(UEs{ueIdx}.MACEntity, channelQualityInfoDL);
    
    % Setup logical channels
    configureLogicalChannel(gNB, ueIdx, rlcChannelConfigStruct);
    configureLogicalChannel(UEs{ueIdx}, ueIdx, rlcChannelConfigStruct);
    
    % Add traffic applications
    ulApp = networkTrafficOnOff('GeneratePacket', true, ...
        'OnTime', simParameters.NumFramesSim/100, 'OffTime', 0, ...
        'DataRate', simParameters.ulAppDataRate(ueIdx));
    UEs{ueIdx}.addApplication(ueIdx, simParameters.LCHConfig.LCID, ulApp);
    
    dlApp = networkTrafficOnOff('GeneratePacket', true, ...
        'OnTime', simParameters.NumFramesSim/100, 'OffTime', 0, ...
        'DataRate', simParameters.dlAppDataRate(ueIdx));
    gNB.addApplication(ueIdx, simParameters.LCHConfig.LCID, dlApp);
    
    fprintf('UE %d created\n', ueIdx);
end

fprintf('\n');

%% ========================================================================
% SECTION 5: CREATE NETWORK SIMULATOR
% ========================================================================

fprintf('CREATING NETWORK SIMULATOR...\n\n');

nrNodes = [{gNB}; UEs];
networkSimulator = hWirelessNetworkSimulator(nrNodes);

fprintf('Network simulator created\n\n');

%% ========================================================================
% SECTION 6: LOGGING
% ========================================================================

fprintf('SETTING UP LOGGING...\n\n');

if ~exist('Results', 'dir')
    mkdir('Results');
end

if enableTraces
    simRLCLogger = hNRRLCLogger(simParameters, lchInfo, networkSimulator, gNB, UEs);
    simSchedulingLogger = hNRSchedulingLogger(simParameters, networkSimulator, gNB, UEs);
end

simParameters.MetricsStepSize = ceil((simParameters.NumFramesSim * 10 * simParameters.SCS / 15) / simParameters.NumMetricsSteps);

ViewPlots = true;
metricsVisualizer = hNRMetricsVisualizer(simParameters, ...
    'EnableSchedulerMetricsPlots', ViewPlots, ...
    'EnableRLCMetricsPlots', ViewPlots, ...
    'LCHInfo', lchInfo, ...
    'NetworkSimulator', networkSimulator, ...
    'GNB', gNB, 'UEs', UEs);

fprintf('Logging configured\n\n');

%% ========================================================================
% SECTION 7: REGISTER RADAR CHANNEL MODEL
% ========================================================================

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║  REGISTERING CUSTOM CHANNEL WITH RADAR INJECTION          ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('Channel model: hAddImpairmentWithRadar\n');
fprintf('This function is called for EVERY packet transmission:\n');
fprintf('  - DL packets: gNB → UE (through CDL-C + RADAR + Noise)\n');
fprintf('  - UL packets: UE → gNB (through CDL-C + RADAR + Noise)\n\n');

fprintf('For each transmission:\n');
fprintf('  1. CDL-C multipath fading applied\n');
fprintf('  2. Radar waveform ADDED to signal\n');
fprintf('  3. AWGN noise added\n');
fprintf('  4. Degraded signal returned to PHY layer\n');
fprintf('  5. PHY calculates SINR from degraded signal\n');
fprintf('  6. CQI reported to scheduler\n');
fprintf('  7. HARQ decides if retransmission needed\n\n');

% THIS IS THE KEY LINE - REGISTER THE CUSTOM CHANNEL FUNCTION
% The function hAddImpairmentWithRadar will be called by the network
% simulator for EVERY transmission (UL and DL)
fprintf('Registering channel model with network simulator...\n\n');

% Create anonymous function that captures radar parameters
radarWaveformHandle = simParameters.radar.waveform;
radarFreqOffsetHandle = simParameters.radar.FreqOffset;
sampleRateHandle = simParameters.radar.SampleRate;

% Register the custom channel function
addChannelModel(networkSimulator, ...
    @(txWaveform) hAddImpairmentWithRadar(txWaveform, ...
        radarWaveformHandle, radarFreqOffsetHandle, sampleRateHandle));



%% ========================================================================
% SECTION 8: RUN SIMULATION
% ========================================================================


simulationTime = simParameters.NumFramesSim * 1e-2;

fprintf('Starting simulation for %.1f seconds...\n\n', simulationTime);

run(networkSimulator, simulationTime);

fprintf('\n Simulation completed!\n\n');



displayPerformanceIndicators(metricsVisualizer);

%% ========================================================================
% SECTION 9: SAVE LOGS
% ========================================================================

fprintf('SAVING SIMULATION LOGS...\n\n');

% Create simulation logs from loggers
simulationLogs = cell(1,1);
logInfo = struct('TimeStepLogs', [], 'SchedulingAssignmentLogs', [], 'RLCLogs', []);

if enableTraces
    try
        [logInfo.TimeStepLogs] = getSchedulingLogs(simSchedulingLogger);
        logInfo.SchedulingAssignmentLogs = getGrantLogs(simSchedulingLogger);
        logInfo.RLCLogs = getRLCLogs(simRLCLogger);
        simulationLogs{1} = logInfo;
        fprintf(' Logs extracted from loggers\n\n');
    catch ME
        fprintf(' Warning: Could not extract logs: %s\n', ME.message);
        fprintf('  Retransmission analysis will be skipped\n\n');
        simulationLogs = {};
    end
else
    fprintf('Tracing disabled - no logs available\n\n');
end

%% ========================================================================
% SECTION 10: RETRANSMISSION ANALYSIS
% ========================================================================

fprintf('ANALYZING RETRANSMISSIONS (HARQ IMPACT)...\n\n');

if ~isempty(simulationLogs) && ~isempty(simulationLogs{1})
try
    simTable = simulationLogs{1}.SchedulingAssignmentLogs();
    numReTxDL = zeros(simParameters.NumUEs, 1);
    numReTxUL = zeros(simParameters.NumUEs, 1);
    numNewTxDL = zeros(simParameters.NumUEs, 1);
    numNewTxUL = zeros(simParameters.NumUEs, 1);
    
    [rows, cols] = size(simTable);
    
    for i = 2:rows-1
        try
            rnti = simTable{i, 1};
            txType = simTable{i, 13};
            direction = simTable{i, 4};
            
            if string(txType) == 'reTx' && string(direction) == 'DL'
                numReTxDL(rnti) = numReTxDL(rnti) + 1;
            elseif string(txType) == 'newTx' && string(direction) == 'DL'
                numNewTxDL(rnti) = numNewTxDL(rnti) + 1;
            end
            
            if string(txType) == 'reTx' && string(direction) == 'UL'
                numReTxUL(rnti) = numReTxUL(rnti) + 1;
            elseif string(txType) == 'newTx' && string(direction) == 'UL'
                numNewTxUL(rnti) = numNewTxUL(rnti) + 1;
            end
        catch
        end
    end
    
    fprintf('DOWNLINK:\n');
    for i = 1:simParameters.NumUEs
        total_dl = numNewTxDL(i) + numReTxDL(i);
        if total_dl > 0
            rate = 100 * numReTxDL(i) / total_dl;
            fprintf('  UE %d: %d new + %d retx = %.1f%% retx rate\n', ...
                i, numNewTxDL(i), numReTxDL(i), rate);
        end
    end
    
    fprintf('\nUPLINK:\n');
    for i = 1:simParameters.NumUEs
        total_ul = numNewTxUL(i) + numReTxUL(i);
        if total_ul > 0
            rate = 100 * numReTxUL(i) / total_ul;
            fprintf('  UE %d: %d new + %d retx = %.1f%% retx rate\n', ...
                i, numNewTxUL(i), numReTxUL(i), rate);
        end
    end
    
    totalNewDL = sum(numNewTxDL);
    totalReTxDL = sum(numReTxDL);
    totalNewUL = sum(numNewTxUL);
    totalReTxUL = sum(numReTxUL);
    
    fprintf('\nTOTAL:\n');
    if (totalNewDL + totalReTxDL) > 0
        dlRate = 100 * totalReTxDL / (totalNewDL + totalReTxDL);
        fprintf('  DL: %.1f%% retransmission rate\n', dlRate);
    end
    if (totalNewUL + totalReTxUL) > 0
        ulRate = 100 * totalReTxUL / (totalNewUL + totalReTxUL);
        fprintf('  UL: %.1f%% retransmission rate\n', ulRate);
    end
    
catch ME
    fprintf('(Could not extract retransmission stats: %s)\n', ME.message);
end
else
    fprintf('No logs available - skipping retransmission analysis\n');
    fprintf('(Enable enableTraces = true for detailed statistics)\n\n');
end



