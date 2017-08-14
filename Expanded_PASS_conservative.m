% Testing the Expanded Parametric Adaptive Spectrum Sensing (PASS) algorithm for
%  dynamic spectrum sensing.
% Expanded algorithm includes transmit scheduling as well as scan
%  scheduling.
% From 2007 paper by D. Datla et al. and 2009 paper by Datla et al.
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
addonMax = 10;                        % Maximum length of a single transmit period
occupied2 = zeros(channels, 1);
vacant2 = zeros(channels, 1);
n1 = zeros(channels, 1);        % Scan period backoff
n2 = zeros(channels, 1);        % Transmit period modifier
samples = zeros(channels, 1);   % # of times each channel sampled by PASS
sweepInst = zeros(length, 1);
transmit = zeros(channels, length);
interfere = zeros(channels, length);

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
for j = 1:length
    for i = 1:channels
        current = M(i, j);
        if A(i, j) == 1            
            samples(i) = samples(i) + 1;
            sweepInst(j) = sweepInst(j) + 1;
            if current == 1
                occupied2(i) = occupied2(i) + 1;
                n1(i) = n1(i) + 1;
                if n1(i) > backoffMax
                    n1(i) = backoffMax;
                end
                n2(i) = 0;
                temp1 = j + 1;
                temp2 = j + n1(i);
                if temp1 > length
                    temp1 = length;
                end
                if temp2 > length
                    temp2 = length;
                end
                A(i, temp1: temp2) = 0;
            elseif current == 0
                vacant2(i) = vacant2(i) + 1;
                n2(i) = n2(i) + 1;
                if n2(i) > addonMax
                    n2(i) = addonMax;
                end
                n1(i) = 0;
                temp1 = j + 1;
                temp2 = j + 1;              % conservative approach
%                 temp2 = j + n2(i);          % non-conservative approach
                if temp1 > length
                    temp1 = length;
                end
                if temp2 > length
                    temp2 = length;
                end
                A(i, temp1: temp2) = 2;
            end 
        elseif A(i, j) == 2
            if current == 1
                interfere(i, j) = 1;
            elseif current == 0
                transmit(i, j) = 1;
            end
        elseif A(i, j) == 0
%             n(i) = 1;
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
spectrumUtil = 100.*(sum(transmit, 2) ./ vacant);
for i = 1:channels
    if vacant(i) == 0
       spectrumUtil(i, :) = 0;
    end
end
spectrumUtilMean = mean(spectrumUtil);
interfereRate = 100.*(sum(interfere, 2)./length);
interfereMean = mean(interfereRate);

M2 = abs(M' - 1);
% figure
% image(M2, 'CDataMapping', 'scaled')
% colormap gray 