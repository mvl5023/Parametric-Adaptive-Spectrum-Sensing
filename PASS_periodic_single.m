% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% Single channel testing, periodic occupancy
% Sweeping across max scan period and occupancy percentage
% From 2007 paper by D. Datla et al.
%-----------------------------------------------------------------------

% Simulation parameters
length = 100000;          % # of samples in occupancy matrix

% Variables for PASS algorithm
%A = ones( channels , t );       % Matrix of time-frequency assignments
backoffMax = 10;                 % Maximum backoff of PASS algorithm
startP = 1;
stopP = 10;
sweepsP = 10;
startQ = 10;
stopQ = 200;
sweepsQ = 20;
reductionTot = zeros(sweepsP, sweepsQ);    % p x q
lossTot = zeros(sweepsP, sweepsQ);
samplesTot = zeros(sweepsP, sweepsQ);
vacanciesTot = zeros(sweepsP, sweepsQ);
vacanciesRatio = zeros(sweepsP, sweepsQ);

for p = linspace(startP, stopP, sweepsP)           % channel occupancy variation
    x = round((p - startP) + 1);
     
    M = [ zeros(1, stopP - p) , ones(1, p) ];
    M = repmat(M, 1, length/10);
    
    % Calculate number of occupied and vacant samples per channel
    occupied = sum(M);
    vacant = length - occupied;
    for q = linspace(startQ, stopQ, sweepsQ)        % length of time-freq assignment matrix row
        y = (q - startQ)/10 + 1;
        A = ones( 1 , q );       % Matrix of time-frequency assignments
        sweeps = round(length / q);             % # of times PASS algorithm runs
        occupied2 = 0;
        vacant2 = 0;
        n = 1;          % Scan period multiplier array
        samples = 0;   % # of times each channel sampled by PASS
        for k = 1:sweeps
            %A = ones( 1 , q );       % Matrix of time-frequency assignments
            for j = 1:q                
                current = (k - 1)*10 + j;
                if A(j) == 1
                    temp = M(current);
                    samples = samples + 1;
                    if temp == 1
                        occupied2 = occupied2 + 1;
                        n = n + 1;
%                         if n > min([backoffMax, q])
%                             n = min([backoffMax, q]);
%                         end
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
                        A(j: j+1) = 1;      % Change from paper's algorithm
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
        vacanciesTot(x, y) = vacant2;
    end  
    
    if p < 10
        vacanciesRatio(x, :) = vacanciesTot(x, :)./((10 - p)*10000); 
    elseif p == 10
        vacanciesRatio(x, :) = 0;
    end
    
end