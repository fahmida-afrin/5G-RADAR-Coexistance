%% NR TDD Symbol-Based Scheduling — Segmented & Commented Setup
% This script configures and runs a small NR (5G) network simulation using
% Communications Toolbox Wireless Network Simulation Library + 5G Toolbox.

%% 0) Environment & Support Package Check
wirelessnetworkSupportPackageCheck;        % Verify WNSL support package
rng('default');                            % Reproducible randomness

%% 1) Primary Scenario Parameters 
simParameters = struct();

% — Simulation horizon —
simParameters.NumFramesSim   = 10;   % # of 10 ms frames (10 -> 0.1 s)

% — Scheduler & duplex mode —
simParameters.SchedulingType = 1;    % 0 = slot-based, 1 = symbol-based
simParameters.DuplexMode     = 1;    % 0 = FDD, 1 = TDD
simParameters.SchedulerStrategy = 'RR'; % 'PF' | 'RR' | 'BestCQI' | 'StaticMCS'

% — Carrier numerology —
simParameters.SCS            = 15;          % kHz (15, 30, 60, ...)
simParameters.NumRBs         = 20;          % Resource Blocks (bandwidth grid)
simParameters.DLBandwidth    = 5e6;         % Hz (informational)
simParameters.ULBandwidth    = 5e6;         % Hz (informational)
simParameters.DLCarrierFreq  = 3.595e9;     % Hz
simParameters.ULCarrierFreq  = 3.595e9;     % Hz
simParameters.NCellID        = 1;           % Physical cell ID

% — Topology —
simParameters.Position       = [0 0 0];     % gNB position (x,y,z) meters
simParameters.NumUEs         = 4;
simParameters.UEPosition     = [            % UE positions (x,y,z) meters
    100  0   0
    150  0   0
    300  0   0
    400  0   0
];

% — TDD DL/UL pattern (Type B mapping for symbol-based by default) —
simParameters.DLULPeriodicity = 5;   % ms per DL/UL period
simParameters.NumDLSlots      = 2;   % full-DL slots at start of each period
simParameters.NumDLSyms       = 8;   % DL symbols at start of transition slot
simParameters.NumULSyms       = 4;   % UL symbols at end of transition slot
simParameters.NumULSlots      = 2;   % full-UL slots at end of each period

% — Symbol-level TTI granularity (only used when SchedulingType = 1) —
simParameters.TTIGranularity   = 4;  % symbols per TTI (e.g., 2/4/7/14)

% — RB allocation limits for NEW transmissions (re-Tx not limited) —
simParameters.RBAllocationLimitUL = 15;  % PUSCH
simParameters.RBAllocationLimitDL = 15;  % PDSCH

% — Reporting/PHY prep —
simParameters.BSRPeriodicity   = 1;      % ms, UL BSR periodicity
simParameters.PUSCHPrepTime    = 200;    % us, PUSCH assignment lead time

% — Synthetic channel quality evolution (for passthrough PHY / link-abstr) —
simParameters.ChannelUpdatePeriodicity = 0.2;  % seconds
simParameters.CQIDelta                 = 1;    % CQI step change
% [distance(m)  maxUL_CQI]
simParameters.CQIvsDistance = [
    200    15
    300    12
    500    10
    1000    8
    1200    7
];

% — DM-RS configuration —
simParameters.DMRSTypeAPosition = 2; % {2,3}
% PUSCH
simParameters.PUSCHDMRSAdditionalPosTypeB = 0;
simParameters.PUSCHDMRSAdditionalPosTypeA = 0;
simParameters.PUSCHDMRSConfigurationType  = 1;
% PDSCH
simParameters.PDSCHDMRSAdditionalPosTypeB = 0;
simParameters.PDSCHDMRSAdditionalPosTypeA = 0;
simParameters.PDSCHDMRSConfigurationType  = 1;

% — Application traffic (per-UE, On-Off always On here) —
% Units are bits per second (bps). Set both UL and DL.
simParameters.dlAppDataRate = [16e5, 30e4, 100e4, 10e4];  % 1.6 Mb/s, 0.3 Mb/s, 1.0 Mb/s, 0.1 Mb/s
simParameters.ulAppDataRate = [16e5, 30e4, 100e4, 10e4];

% — Logging & visualization —
simParameters.CQIVisualization = false;   % heatmap of CQI per RB
simParameters.RBVisualization  = false;   % RB assignment grid
simParameters.NumMetricsSteps  = 20;      % UI update steps across sim
enableTraces                   = true;    % log RLC/MAC traces to MATs

% — Output file names —
parametersLogFile  = 'simParameters';
simulationLogFile  = 'simulationLogs';
simulationMetricsFile = 'simulationMetrics';

%% 2) Derived Parameters & Sanity Checks
% Slots per 10 ms frame = 10 * (SCS/15 kHz)
numSlotsSim = (simParameters.NumFramesSim * 10 * simParameters.SCS) / 15;

% Mapping type based on SchedulingType
if simParameters.SchedulingType == 1
    simParameters.PUSCHMappingType = 'B';
    simParameters.PDSCHMappingType = 'B';
else
    simParameters.PUSCHMappingType = 'A';
    simParameters.PDSCHMappingType = 'A';
end

% Metrics step guardrail
simParameters.MetricsStepSize = max(1, ceil(numSlotsSim / simParameters.NumMetricsSteps));

% RLC & Logical channel config
numLogicalChannels = 1; 
simParameters.LCHConfig.LCID = 4;
simParameters.RLCConfig.EntityType = 2;  % 0:UM DL, 1:UM UL, 2:UM bidi, 3:AM
simParameters.maxRLCSDULength  = 9000;   % per TS 38.323

% Build RLC logger info per UE
validateattributes(simParameters.UEPosition, {'numeric'}, ...
    {'nonempty','real','nrows',simParameters.NumUEs,'ncols',3,'finite'}, ...
    'simParameters.UEPosition','UEPosition');

lchInfo = repmat(struct('RNTI',[], 'LCID',[], 'EntityDir',[]), [simParameters.NumUEs 1]);
for idx = 1:simParameters.NumUEs
    lchInfo(idx).RNTI      = idx;
    lchInfo(idx).LCID      = simParameters.LCHConfig.LCID;
    lchInfo(idx).EntityDir = simParameters.RLCConfig.EntityType;
end

% Logical Channel config struct used by both gNB & UEs
rlcChannelConfigStruct = struct( ...
    'LCGID', 1, ...      % logical channel group
    'Priority', 1, ...   % higher = more priority
    'PBR', 8, ...        % kB/s
    'BSD', 10, ...       % ms
    'EntityType', simParameters.RLCConfig.EntityType, ...
    'LogicalChannelID', simParameters.LCHConfig.LCID);

% Validate app data rates
validateattributes(simParameters.dlAppDataRate, {'numeric'}, ...
    {'nonempty','vector','numel',simParameters.NumUEs,'finite','>',0}, 'dlAppDataRate','dlAppDataRate');
validateattributes(simParameters.ulAppDataRate, {'numeric'}, ...
    {'nonempty','vector','numel',simParameters.NumUEs,'finite','>',0}, 'ulAppDataRate','ulAppDataRate');

%% 3) CQI Initialization (UL & DL) with Euclidean Distance Mapping
% Compute max UL CQI per UE from UE–gNB distance using the mapping.
maxUECQIs = zeros(simParameters.NumUEs, 1);
for ueIdx = 1:simParameters.NumUEs
    d = norm(simParameters.UEPosition(ueIdx,:) - simParameters.Position); % meters
    % Find the first mapping row whose distance threshold exceeds UE distance
    row = find(simParameters.CQIvsDistance(:,1) >= d, 1, 'first');
    if isempty(row)
        row = size(simParameters.CQIvsDistance,1); % farthest bucket
    end
    maxUECQIs(ueIdx) = simParameters.CQIvsDistance(row, 2);
end

% Random per-RB CQI subject to per-UE maxima (UL = DL initially)
simParameters.InitialChannelQualityUL = zeros(simParameters.NumUEs, simParameters.NumRBs);
simParameters.InitialChannelQualityDL = zeros(simParameters.NumUEs, simParameters.NumRBs);
for ueIdx = 1:simParameters.NumUEs
    cqiUL = randi([1, maxUECQIs(ueIdx)], 1, simParameters.NumRBs);
    simParameters.InitialChannelQualityUL(ueIdx,:) = cqiUL;
    simParameters.InitialChannelQualityDL(ueIdx,:) = cqiUL; % start equal
end

%% 4) Node Creation: gNB + UEs, PHY/MAC wiring, Logical Channels & Apps
% — gNB —
gNB = hNRGNB(simParameters);
switch simParameters.SchedulerStrategy
    case 'RR'
        scheduler = hNRSchedulerRoundRobin(simParameters);
    case 'PF'
        scheduler = hNRSchedulerProportionalFair(simParameters);
    case 'BestCQI'
        scheduler = hNRSchedulerBestCQI(simParameters);
    case 'StaticMCS'
        scheduler = hNRSchedulerStaticMCS(simParameters);
    otherwise
        error('Unsupported SchedulerStrategy: %s', simParameters.SchedulerStrategy);
end
addScheduler(gNB, scheduler);

% PHY selection (comment/uncomment as needed)
% gNB.PhyEntity = hNRGNBPassThroughPhy(simParameters); % Link-abstraction only
gNB.PhyEntity = hNRGNBPhy(simParameters);              % Full PHY
configurePhy(gNB, simParameters);
setPhyInterface(gNB);

% — UEs —
UEs = cell(simParameters.NumUEs, 1);
for ueIdx = 1:simParameters.NumUEs
    simParameters.Position = simParameters.UEPosition(ueIdx,:); % temp: per-UE
    UEs{ueIdx} = hNRUE(simParameters, ueIdx);

    % Uplink/Downlink PHY for UE (select one)
    % UEs{ueIdx}.PhyEntity = hNRUEPassThroughPhy(simParameters, ueIdx);
    UEs{ueIdx}.PhyEntity   = hNRUEPhy(simParameters, ueIdx);
    configurePhy(UEs{ueIdx}, simParameters);
    setPhyInterface(UEs{ueIdx});

    % Seed scheduler with initial CQIs (UL & DL)
    updateChannelQualityUL(gNB.MACEntity.Scheduler, struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityUL(ueIdx,:)));
    updateChannelQualityDL(gNB.MACEntity.Scheduler, struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityDL(ueIdx,:)));

    % For UE-side packet error probability estimation (DL)
    updateChannelQualityDL(UEs{ueIdx}.MACEntity, struct('RNTI', ueIdx, 'CQI', simParameters.InitialChannelQualityDL(ueIdx,:)));

    % Logical channels on both ends
    configureLogicalChannel(gNB, ueIdx, rlcChannelConfigStruct);
    configureLogicalChannel(UEs{ueIdx}, ueIdx, rlcChannelConfigStruct);

    % Application traffic (On-Off with always-on period across sim horizon)
    onTimeSec = simParameters.NumFramesSim * 10e-3; % frames * 10ms
    ulApp = networkTrafficOnOff('GeneratePacket', true, 'OnTime', onTimeSec, 'OffTime', 0, 'DataRate', simParameters.ulAppDataRate(ueIdx));
    UEs{ueIdx}.addApplication(ueIdx, simParameters.LCHConfig.LCID, ulApp);

    dlApp = networkTrafficOnOff('GeneratePacket', true, 'OnTime', onTimeSec, 'OffTime', 0, 'DataRate', simParameters.dlAppDataRate(ueIdx));
    gNB.addApplication(ueIdx, simParameters.LCHConfig.LCID, dlApp);
end

%% 5) Simulator, Loggers, & Visualizers
nrNodes = [{gNB}; UEs];
networkSimulator = hWirelessNetworkSimulator(nrNodes);

if enableTraces
    % RLC metrics logged per-slot
    simRLCLogger        = hNRRLCLogger(simParameters, lchInfo, networkSimulator, gNB, UEs);
    simSchedulingLogger = hNRSchedulingLogger(simParameters, networkSimulator, gNB, UEs);

    % Optional CQI/RB grid visualizer
    if simParameters.CQIVisualization || simParameters.RBVisualization
        gridVisualizer = hNRGridVisualizer(simParameters, 'MACLogger', simSchedulingLogger); 
    end
end

% Achieved rate / spectral efficiency plots
ViewPlots = true;
metricsVisualizer = hNRMetricsVisualizer(simParameters, ...
    'EnableSchedulerMetricsPlots', ViewPlots, ...
    'EnableRLCMetricsPlots', ViewPlots, ...
    'LCHInfo', lchInfo, ...
    'NetworkSimulator', networkSimulator, ...
    'GNB', gNB, 'UEs', UEs);

%% 6) Channel/Impairments Hook (Optional)
% Plug in a custom impairment/channel function if available
% Signature expected by addChannelModel(): @(ns, ...) hAddImpairment(ns, ...)
try
    addChannelModel(networkSimulator, @hAddImpairment);
catch ME
    warning('addChannelModel skipped (%s). Proceeding without custom impairment.', ME.message);
end

%% 7) Run Simulation
simulationTime = simParameters.NumFramesSim * 10e-3; % seconds
run(networkSimulator, simulationTime);

displayPerformanceIndicators(metricsVisualizer);
metrics = getMetrics(metricsVisualizer);
save(simulationMetricsFile, 'metrics');

%% 8) Persist Traces & Parameters (If Enabled)
if enableTraces
    fprintf('Saving logs...\n');

    simulationLogs = cell(1,1);
    logInfo = struct('TimeStepLogs',[], 'SchedulingAssignmentLogs',[], 'RLCLogs',[]);
    logInfo.TimeStepLogs             = getSchedulingLogs(simSchedulingLogger);
    logInfo.SchedulingAssignmentLogs = getGrantLogs(simSchedulingLogger);
    logInfo.RLCLogs                  = getRLCLogs(simRLCLogger);
    simulationLogs{1} = logInfo;

    save(simulationLogFile, 'simulationLogs');
    save(parametersLogFile, 'simParameters');

    dt = datestr(now,'yymmdd-HHMMSS');
    newFolderName = fullfile('Results', dt, filesep);
    if ~exist(newFolderName, 'dir'); mkdir(newFolderName); end

    save(fullfile(newFolderName, sprintf('%s_%s.mat', dt, simulationMetricsFile)), 'simulationLogs');
    save(fullfile(newFolderName, sprintf('%s_%s.mat', dt, parametersLogFile)), 'simParameters');

    fprintf('Logs saved under: %s\n', newFolderName);
end


