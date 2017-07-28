% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%-----------------------------------------------------------------------

% Simulation parameters
channels = 10;          % # of channels in spectrum band
length = 10000;          % # of samples in occupancy matrix
m = 0.6;                % constant coefficient for exponential distr.
b = 0.02;               % offset for exponential distr.
t = 10;                 % Spectrum sweep max period

% Generate test matrix of randomly generated spectrum occupancy data using
% exponential distribution
M = spectrum_occ_exp( channels, length, m, b);

% Calculate number of occupied and vacant samples per channel
occupied = sum(M, 2);
vacant = length - occupied;

% Variables for PASS algorithm
%A = ones( channels , t );       % Matrix of time-frequency assignments
occupied2 = zeros(channels, 1);
vacant2 = zeros(channels, 1);
n = ones(channels, 1);          % Scan period multiplier array
samples = zeros(channels, 1);   % # of times each channel sampled by PASS
sweeps = length / t;            % # of times PASS algorithm runs
backoffMax = 10;                % Maximum backoff of PASS algorithm

% Run PASS algorithm
for k = 1:sweeps
    A = ones( channels , t );       % Matrix of time-frequency assignments
    for j = 1:t
        for i = 1:channels
            current = (k - 1)*10 + j;
            if A(i, j) == 1
                temp = M(i, current);
                samples(i) = samples(i) + 1;
                if temp == 1
                    occupied2(i) = occupied2(i) + 1;
                    n(i) = n(i) + 1;
                    if n(i) > backoffMax
                        n(i) = backoffMax;
                    end
                    temp2 = j + n(i) - 1;
                    if temp2 > t
                        temp2 = t;
                    end
                    A(i, j: temp2) = 0;
                elseif temp == 0
                    vacant2(i) = vacant2(i) + 1;
                    n(i) = 1;
                end 
            end
        end
    end
end

% Calculate metrics
reduction = samples ./ length;       % Scans as fraction of total length
reductionTot = sum(samples) / (channels * length);       
loss = vacant2 ./ vacant;
lossTot = sum(vacant2)/sum(vacant);

M = 50.* M;
figure
image(M)
colormap hot
   








