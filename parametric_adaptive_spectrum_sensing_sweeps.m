% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%  * Multi-channel testing, exponential distribution
%  * Sweeping across max scan period and exponential distr. coefficient
%-----------------------------------------------------------------------

% Simulation parameters
channels = 10;          % # of channels in spectrum band
length = 100000;          % # of samples in occupancy matrix
b = 0.02;               % offset for exponential distr.

% Variables for PASS algorithm
A = ones( channels , length );       % Matrix of time-frequency assignments
startP = 0.5;
stopP = 1.2;
sweepsP = 8;
startQ = 10;
stopQ = 200;
sweepsQ = 20;
n1_def = 1;
n2_def = 1;
reductionTot = zeros(sweepsP, sweepsQ, channels);
samplesTot = zeros(sweepsP, sweepsQ, channels);
vacanciesTot = zeros(sweepsP, sweepsQ, channels);
vacanciesTot2 = zeros(sweepsP, sweepsQ, channels);

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
for p = linspace(startP, stopP, sweepsP)           % exponential distribution coefficient
    x = round((p - startP)/0.1 + 1);
    m = p;                % constant coefficient for exponential distr.
    
    % Generate test matrix of randomly generated spectrum occupancy data using
    % exponential distribution
    M = spectrum_occ_exp( channels, length, m, b);

    % Calculate number of occupied and vacant samples per channel
    occupied = sum(M, 2);
    vacant = length - occupied;
    for q = linspace(startQ, stopQ, sweepsQ)        % length of time-freq assignment matrix row
        y = (q - startQ)/10 + 1;
        occupied2 = zeros(channels, 1);
        vacant2 = zeros(channels, 1);
        n1 = n1_def .* ones(channels, 1);          % Number of scan periods removed from A
        n2 = n2_def .* ones(channels, 1);          % Number of scan periods added to A
        samples = zeros(channels, 1);   % # of times each channel sampled by PASS
        for j = 1:length         % Sweep through columns of A (samples)
            for i = 1:channels    % Sweep through rows of A (channels)
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
                        n2(i) = n2_def;
                    elseif temp == 0
                        % Original algorithm--------------------------
                        vacant2(i) = vacant2(i) + 1;
                        n1(i) = n1_def;
                        %---------------------------------------------
                        %*********************************************
                        % Modified with restoration-------------------
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
                elseif A(i, j) == 0
                    n1(i) = 1;
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