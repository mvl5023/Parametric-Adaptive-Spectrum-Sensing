% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%  * Single channel testing, exponential distribution
%  * Sweeping across max scan period and exponential distr. coefficient
%-----------------------------------------------------------------------

% Simulation parameters
length = 100000;          % # of samples in occupancy matrix
b = 0.02;               % offset for exponential distr.

% Variables for PASS algorithm
%A = ones( channels , t );       % Matrix of time-frequency assignments
%backoffMax = 10;                 % Maximum backoff of PASS algorithm
startP = 0.3;
stopP = 0.6;
sweepsP = 7;
startQ = 10;
stopQ = 200;
sweepsQ = 20;
reductionTot = zeros(sweepsP, sweepsQ);    % p x q
lossTot = zeros(sweepsP, sweepsQ);
samplesTot = zeros(sweepsP, sweepsQ);
vacanciesTot = zeros(sweepsP, sweepsQ);
vacanciesTot2 = zeros(sweepsP, sweepsQ);

for p = linspace(startP, stopP, sweepsP)           % exponential distribution coefficient
    x = round((p - startP)/0.05 + 1);
    m = p;                % constant coefficient for exponential distr.
    
    % Generate test matrix of randomly generated spectrum occupancy data using
    % exponential distribution
    M = spectrum_occ_exp( 1, length, m, b);

    % Calculate number of occupied and vacant samples per channel
    occupied = sum(M);
    vacant = length - occupied;
    for q = linspace(startQ, stopQ, sweepsQ)        % length of time-freq assignment matrix row
        y = (q - startQ)/10 + 1;
        A = ones( 1 , q );       % Matrix of time-frequency assignments
        sweeps = floor(length / q);             % # of times PASS algorithm runs
        lengthB = sweeps*q;                     % Number of samples included
        % Calculate number of occupied and vacant samples per channel
        occupied = sum(M(1:lengthB));
        vacant = lengthB - occupied;
        occupied2 = 0;
        vacant2 = 0;
        n = 1;          % Scan period multiplier array
        samples = 0;   % # of times each channel sampled by PASS
        for k = 1:sweeps
            %A = ones( 1 , q );       % Matrix of time-frequency assignments
            for j = 1:q                
                current = (k - 1)*q + j;
                if A(j) == 1
                    temp = M(current);
                    samples = samples + 1;
                    if temp == 1
                        occupied2 = occupied2 + 1;
                        n = n + 1;
                        if n > q
                            n = q;
                        end
                        temp2 = j + n - 1;
                        if temp2 > q
                            temp2 = q;
                        end
                        A(j: temp2) = 0;
                    elseif temp == 0
                        vacant2 = vacant2 + 1;
                        n = 1;
                    end 
                elseif A(j) == 0
                    n = 1;
                end                
            end
        end

        % Calculate metrics     
        samplesTot(x, y) = samples;
        reductionTot(x, y) = samples / length;       
        vacanciesTot(x, y) = vacant;
        vacanciesTot2(x, y) = vacant2;
    end    
end

vacanciesRatio = vacanciesTot2 ./ vacanciesTot;
vacanciesRatio(sweepsP, :) = 1;
opt = vacanciesRatio ./ reductionTot;