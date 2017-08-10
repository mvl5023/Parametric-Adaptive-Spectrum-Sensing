% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%-----------------------------------------------------------------------

% Simulation parameters
channels = 3300;          % # of channels in spectrum band
length = 1000;          % # of samples in occupancy matrix
m = 1.2;                % constant coefficient for exponential distr.
b = 0.02;               % offset for exponential distr.

% Generate test matrix of randomly generated spectrum occupancy data using
% exponential distribution
M = spectrum_occ( channels , length );
%M = spectrum_occ_exp( channels, length, m, b);
M(75:600, :) = 1;
M(675:1200, :) = 1;
M(1275:1800, :) = 1;
M(2175:2700, :) = 1;
M(2775:3300, :) = 1;

% Calculate number of occupied and vacant samples per channel
occupied = sum(M, 2);
vacant = length - occupied;

% Variables for PASS algorithm
A = ones( channels , length );       % Matrix of time-frequency assignments
backoffMax = 100;                     % Maximum backoff of PASS algorithm
occupied2 = zeros(channels, 1);
vacant2 = zeros(channels, 1);
n = ones(channels, 1);          % Scan period multiplier array
samples = zeros(channels, 1);   % # of times each channel sampled by PASS
sweepInst = zeros(length, 1);

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
for j = 1:length
    for i = 1:channels
        if A(i, j) == 1
            temp = M(i, j);
            samples(i) = samples(i) + 1;
            sweepInst(j) = sweepInst(j) + 1;
            if temp == 1
                occupied2(i) = occupied2(i) + 1;
                n(i) = n(i) + 1;
                if n(i) > backoffMax
                    n(i) = backoffMax;
                end
                temp2 = j + n(i) - 1;
                if temp2 > length
                    temp2 = length;
                end
                A(i, j: temp2) = 0;
            elseif temp == 0
                vacant2(i) = vacant2(i) + 1;
                n(i) = 1;
            end 
        elseif A(i, j) == 0
%              n(i) = 1;
        end
    end
end 

% Calculate metrics
reduction = samples / length;
reductionTot = sum(samples) / (channels * length);       
vacanciesRatio = vacant2 ./ vacant;
vacanciesRatioMean = 100.*(sum(vacant2)/sum(vacant));
opt = vacanciesRatio ./ reduction;
efficiency = 100.*(length - samples)./length;
efficiencyMean = mean(efficiency);
sensingT = 100.*samples / length;
sweepInst = 100.*sweepInst/channels;

% M2 = abs(M' - 1);
% figure
% image(M2, 'CDataMapping', 'scaled')
% colormap gray 