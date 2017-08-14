% Testing the Expanded Parametric Adaptive Spectrum Sensing (PASS) algorithm for
%  dynamic spectrum sensing.
% Expanded algorithm includes transmit scheduling as well as scan
%  scheduling.
% From 2007 paper by D. Datla et al. and 2009 paper by Datla et al.
%  * Single channel testing, periodic occupancy
%  * Sweeping across max backoff and occupancy percentage
%-----------------------------------------------------------------------

% Simulation parameters
length = 200000;          % # of samples in occupancy matrix

% Variables for PASS algorithm
A = ones( 1 , length );       % Matrix of time-frequency assignments
a = 2;                     % base for exponential backoff
startP = 0;
stopP = 100;
sweepsP = 101;
startQ = 10;
stopQ = 200;
sweepsQ = 20;
n1_def = 0;
n2_def = 0;
reductionTot = zeros(sweepsP, sweepsQ);    % p x q
lossTot = zeros(sweepsP, sweepsQ);
samplesTot = zeros(sweepsP, sweepsQ);
vacanciesTot = zeros(sweepsP, sweepsQ);
vacanciesTot2 = zeros(sweepsP, sweepsQ);
occupTot = zeros(sweepsP, sweepsQ);
occupTot2 = zeros(sweepsP, sweepsQ);
transmitTot = zeros(sweepsP, sweepsQ);
interfereTot = zeros(sweepsP, sweepsQ);

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
% Sweep channel occupancy percent
for p = linspace(startP, stopP, sweepsP)           
    x = p + 1;
    
%     % Generate test array of single channel with periodic occupancy
%     M = [ zeros(1, stopP - p) , ones(1, p) ];
%     M = repmat(M, 1, length/100);

    % Generate test array of single channel with random occupancy
    M = zeros(1, length);
    for k = 1:length
        roll = stopP * rand;
        if p >= roll
            M(k) = 1;
        elseif p < roll
            M(k) = 0;
        end
    end

    % Calculate number of occupied and vacant samples per channel
    occupied = sum(M);
    vacant = length - occupied;
    
    % Sweep maximum backoff
    for q = linspace(startQ, stopQ, sweepsQ)      
        y = (q - startQ)/10 + 1;   
        occupied2 = 0;
        vacant2 = 0;
        transmit = 0;
        interfere = 0;
        
        % Initiate run-specific variables
        A = ones( 1 , length );       % Matrix of time-frequency assignments
        n1 = n1_def;            % Number of scan periods removed from A
        n2 = n2_def;            % Number of scan periods added to A
        samples = 0;            % Number of times each channel sampled by PASS
        
        % Sweep samples
        for j = 1:length  
            current = M(j);
            if A(j) == 1                
                samples = samples + 1;
                if current == 1
                    occupied2 = occupied2 + 1;
                    n1 = n1 + 1;
                    % Max backoff set by user parameter
                    if n1 > q
                        n1 = q;
                    end
                    n2 = n2_def;
                    temp1 = j + 1;
                    % Backoff 1 option
                    %----------------------------------------------------------------
%                     temp2 = j + n1;                 % linear backoff
                    temp2 = j + power(a, n1);    % exponential backoff
                    %----------------------------------------------------------------
                    if temp1 > length
                        temp1 = length;
                    end
                    if temp2 > length
                        temp2 = length;
                    end                        
                    A(temp1: temp2) = 0;
                elseif current == 0
                    vacant2 = vacant2 + 1;
                    n2 = n2 + 1;
                    if n2 > q
                        n2 = q;
                    end
                    n1 = n1_def;
                    temp1 = j + 1;
                    % Backoff 2 option
                    %--------------------------------------------------------------------
%                     temp2 = j + 1;                % conservative approach
%                     temp2 = j + n2;             % linear non-conservative approach
                    temp2 = j + power(a, n2);     % exponential non-conservative approach
                    %--------------------------------------------------------------------
                    if temp1 > length
                        temp1 = length;
                    end
                    if temp2 > length
                        temp2 = length;
                    end
                    A(temp1: temp2) = 2;  
                end
            elseif A(j) == 2
                if current == 1
                    interfere = interfere + 1;
                elseif current == 0
                    transmit = transmit + 1; 
                end               
            elseif A(j) == 0
                % Do nothing
            end                
        end

        % Calculate metrics     
        samplesTot(x, y) = samples;
        reductionTot(x, y) = 100*samples / length;  
        vacanciesTot(x, y) = vacant;
        vacanciesTot2(x, y) = vacant2;
        occupTot(x, y) = occupied;
        occupTot2(x, y) = occupied2;
        transmitTot(x, y) = transmit;
        interfereTot(x, y) = interfere;
    end   
end

% Data processing
efficiency = 100.*(length - samplesTot)./(length);
efficiencyMean = mean(efficiency, 1);
vacanciesRatio = 100.*vacanciesTot2 ./ vacanciesTot;
vacanciesRatio(sweepsP, :) = 100;
opt = vacanciesRatio ./ reductionTot;
spectrumUtil = 100.*(transmitTot ./ vacanciesTot);
spectrumUtil(sweepsP, :) = 0;
spectrumUtilMean = mean(spectrumUtil, 1);
interfereRate = 100.*(interfereTot./length);
interfereMean = mean(interfereRate, 1);

Ha = [efficiency(:, 1), spectrumUtil(:, 1), interfereRate(:, 1)];