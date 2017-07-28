% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% Sweeping across max scan period and exponential distr. coefficient
% From 2007 paper by D. Datla et al.
%-----------------------------------------------------------------------

% Simulation parameters
channels = 10;          % # of channels in spectrum band
length = 100000;          % # of samples in occupancy matrix
b = 0.02;               % offset for exponential distr.

% Variables for PASS algorithm
%A = ones( channels , t );       % Matrix of time-frequency assignments
%backoffMax = 10;                 % Maximum backoff of PASS algorithm
reductionTot = zeros(17, 18);    % p x q
lossTot = zeros(17, 18);
samplesTot = zeros(17, 18);

for p = 0.4:0.1:2           % exponential distribution coefficient
    x = round((p - 0.4)/0.1 + 1);
    m = p;                % constant coefficient for exponential distr.
    % Generate test matrix of randomly generated spectrum occupancy data using
    % exponential distribution
    M = spectrum_occ_exp( channels, length, m, b);

    % Calculate number of occupied and vacant samples per channel
    occupied = sum(M, 2);
    vacant = length - occupied;
    for q = 10:10:200        % length of time-freq assignment matrix row
        y = (q - 10)/10 + 1;
        sweeps = round(length / q);             % # of times PASS algorithm runs
        occupied2 = zeros(channels, 1);
        vacant2 = zeros(channels, 1);
        n = ones(channels, 1);          % Scan period multiplier array
        samples = zeros(channels, 1);   % # of times each channel sampled by PASS
        for k = 1:sweeps
            A = ones( channels , q );       % Matrix of time-frequency assignments
            for j = 1:q
                for i = 1:channels
                    current = (k - 1)*10 + j;
                    if A(i, j) == 1
                        temp = M(i, current);
                        samples(i) = samples(i) + 1;
                        if temp == 1
                            occupied2(i) = occupied2(i) + 1;
                            n(i) = n(i) + 1;
                            if n(i) > q
                                n(i) = q;
                            end
                            temp2 = j + n(i) - 1;
                            if temp2 > q
                                temp2 = q;
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
        samplesTot(x, y) = sum(samples);
        reductionTot(x, y) = sum(samples) / (channels * length);       
        loss = vacant2 ./ vacant;
        lossTot(x, y) = sum(vacant2)/sum(vacant);
    end    
end

xlswrite('PASS_sweep1_samples', samplesTot)
xlswrite('PASS_sweep1_reduction', reductionTot)
xlswrite('PASS_sweep1_loss', lossTot)