% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%  * Multi-channel testing, exponential distribution
%  * Sweeping across max scan period and max backoff
%-----------------------------------------------------------------------

% Simulation parameters
channels = 10;          % # of channels in spectrum band
length = 100000;          % # of samples in occupancy matrix
m = 0.5;                % constant coefficient for exponential distr.
b = 0.02;               % offset for exponential distr.

% Variables for PASS algorithm
A = ones( channels , length );       % Matrix of time-frequency assignments
startB = 10;                    % Minimum backoff of PASS algorithm  
intervalB = 10;                 % Increment size of max backoff
startQ = 10;                    % Minimum length of A
stopQ = 200;                    % Maximum length of A
intervalQ = 10;                  % Increment size of length of A
sweepsQ = 20;                   % Number of values in sweep
n1_def = 1;
n2_def = 1;
reductionTot = zeros(sweepsQ, sweepsQ, channels);
samplesTot = zeros(sweepsQ, sweepsQ, channels);
vacanciesTot = zeros(sweepsQ, sweepsQ, channels);
vacanciesTot2 = zeros(sweepsQ, sweepsQ, channels);

% Generate test matrix of randomly generated spectrum occupancy data using
% exponential distribution
M = spectrum_occ_exp( channels, length, m, b);

% Calculate number of occupied and vacant samples per channel
occupied = sum(M, 2);
vacant = length - occupied;

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
for q = linspace(startQ, stopQ, sweepsQ)           % Width of A
    x = (q - startQ)/intervalQ + 1;             
    sweeps = floor(length / q);             % # of times PASS algorithm runs
    for b = startB : intervalB : q        % Maximum backoff
        y = (b - startB)/10 + 1;
        occupied2 = zeros(channels, 1);
        vacant2 = zeros(channels, 1);
        %A = ones( channels , q );       % Matrix of time-frequency assignments
        n1 = n1_def .* ones(channels, 1);          % Number of scan periods removed from A
        n2 = n2_def .* ones(channels, 1);          % Number of scan periods added to A
        samples = zeros(channels, 1);   % Number of times each channel sampled by PASS
        for j = 1:length      % Sweep through columns of M (samples)
            for i = 1:channels       % Sweep through rows of A (channels)
                if A(i, j) == 1         % Scan
                    temp = M(i, current);
                    samples(i) = samples(i) + 1;
                    if temp == 1
                        occupied2(i) = occupied2(i) + 1;
                        n1(i) = n1(i) + 1;
                        if n1(i) > b
                            n1(i) = b;
                        end
                        temp2 = j + n1(i) - 1;
                        if temp2 > q
                            temp2 = q;
                        end
                        A(i, j: temp2) = 0;
                        n2(i) = n2_def;
                    elseif temp == 0
                        % Original algorithm
                        vacant2(i) = vacant2(i) + 1;
                        n1(i) = n1_def;
                        % --------------------------------------------
                        %*********************************************
                        % Modified with restoration
%                         vacant2(i) = vacant2(i) + 1;
%                         n2(i) = n2(i) + 1;
%                         if n2(i) > q
%                             n2(i) = q;
%                         end
%                         temp2 = j + n2(i) - 1;
%                         if temp2 > q
%                             temp2 = q;
%                         end
%                         A(i, j: temp2) = 1;      % Change from paper's algorithm
%                         n1(i) = n1_def;
                        %---------------------------------------------
                    end
                elseif A(i, j) == 0     % No Scan
                    n1(i) = n1_def;
                end
            end
        end

        % Calculate metrics
        samplesTot(x, y, :) = samples(:);
        vacanciesTot(x, y, :) = vacant(:);
        vacanciesTot2(x, y, :) = vacant2(:);        
        reductionTot(x, y, :) = samples(:) ./ length;     
    end 
end

samplesSum = sum(samplesTot, 3);
vacanciesSum = sum(vacanciesTot, 3);
vacanciesSum2 = sum(vacanciesTot2, 3);
reductionMean = mean(reductionTot, 3);

vacanciesRatio = vacanciesSum2 ./ vacanciesSum;
opt = vacanciesRatio ./ reductionMean;

% xlswrite('PASS_sweep1_samples', samplesTot)
% xlswrite('PASS_sweep1_reduction', reductionTot)
% xlswrite('PASS_sweep1_loss', lossTot)