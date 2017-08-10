% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%  * Multi-channel testing, exponential distribution
%  * Sweeping across max scan period and max backoff
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
startQ = 10;
stopQ = 100;
sweepsQ = 19;
samples = zeros(channels, 1);   % # of times each channel sampled by PASS
reductionTot = zeros(sweepsQ, channels);
samplesTot = zeros(sweepsQ, channels);
vacanciesTot = zeros(sweepsQ, channels);
vacanciesTot2 = zeros(sweepsQ, channels);
occupTot2 = zeros(sweepsQ, channels);

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
% Sweep maximum backoff
for q = linspace(startQ, stopQ, sweepsQ)       
    x = (q - startQ)/5 + 1;
    occupied2 = zeros(channels, 1);
    vacant2 = zeros(channels, 1);
    
    % Initiate run-specific variables
    A = ones( channels , length );       % Matrix of time-frequency assignments
    n1 = ones(channels, 1);           % Number of scan periods removed from A
    n2 = ones(channels, 1);           % Number of scan periods added to A
    samples = zeros(channels, 1);               % Number of times each channel sampled by PASS
    
    % Sweep samples
    for j = 1:length  
        % Sweep channels
        for i = 1:channels    
            if A(i, j) == 1
                temp = M(i, j);
                samples(i) = samples(i) + 1;
                if temp == 1
                    occupied2(i) = occupied2(i) + 1;
                    n1(i) = n1(i) + 1;
                    % Max backoff set by user parameter
                    if n1(i) > q
                        n1(i) = q;
                    end
                    temp2 = j + n1(i) - 1;
                    if temp2 > length
                        temp2 = length;
                    end
                    A(i, j: temp2) = 0;
                    n2(i) = 1;
                elseif temp == 0
                    vacant2(i) = vacant2(i) + 1;
                    n1(i) = 1;
                end
            end
        end
    end
    % Calculate metrics
    samplesTot(x, :) = samples(:);
    vacanciesTot(x, :) = vacant(:);
    vacanciesTot2(x, :) = vacant2(:); 
    occupTot2(x, :) = occupied2(:);
    reductionTot(x, :) = samples(:) ./ length;      
end 

samplesSum = sum(samplesTot, 2);
vacanciesSum = sum(vacanciesTot, 2);
vacanciesSum2 = sum(vacanciesTot2, 2);
reductionMean = mean(reductionTot, 2);

efficiency = 100.*((length*channels) - samplesSum)./(length*channels);
vacanciesRatio = 100.*(vacanciesSum2 ./ vacanciesSum);
opt = vacanciesRatio ./ (100*reductionMean);

% M2 = abs(M' - 1);
% figure
% image(M2, 'CDataMapping', 'scaled')
% colormap gray 