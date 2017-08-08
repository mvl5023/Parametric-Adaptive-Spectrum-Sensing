% Testing the Parametric Adaptive Spectrum Sensing (PASS) algorithm for
% dynamic spectrum sensing
% From 2007 paper by D. Datla et al.
%  * Single channel testing, periodic occupancy
%  * Sweeping across max backoff and occupancy percentage
%-----------------------------------------------------------------------

% Simulation parameters
length = 100000;          % # of samples in occupancy matrix

% Variables for PASS algorithm
A = ones( 1 , length );       % Matrix of time-frequency assignments
startP = 0;
stopP = 10;
sweepsP = 11;
startQ = 10;
stopQ = 200;
sweepsQ = 20;
n1_def = 1;
n2_def = 1;
reductionTot = zeros(sweepsP, sweepsQ);    % p x q
lossTot = zeros(sweepsP, sweepsQ);
samplesTot = zeros(sweepsP, sweepsQ);
vacanciesTot = zeros(sweepsP, sweepsQ);
vacanciesTot2 = zeros(sweepsP, sweepsQ);
occupTot = zeros(sweepsP, sweepsQ);
occupTot2 = zeros(sweepsP, sweepsQ);

%------------------------------------------------------------------------
% PASS algorithm
%------------------------------------------------------------------------
% Sweep channel occupancy percent
for p = linspace(startP, stopP, sweepsP)           
    x = p + 1;
    % Generate test matrix of single channel with periodic occupancy
    M = [ zeros(1, stopP - p) , ones(1, p) ];
    M = repmat(M, 1, length/10);
    
    % Calculate number of occupied and vacant samples per channel
    occupied = sum(M);
    vacant = length - occupied;
    
    % Sweep maximum backoff
    for q = linspace(startQ, stopQ, sweepsQ)      
        y = (q - startQ)/10 + 1;   
        occupied2 = 0;
        vacant2 = 0;
        
        % Initiate run-specific variables
        A = ones( 1 , length );       % Matrix of time-frequency assignments
        n1 = n1_def;            % Number of scan periods removed from A
        n2 = n2_def;            % Number of scan periods added to A
        samples = 0;            % Number of times each channel sampled by PASS
        
        % Sweep samples
        for j = 1:length              
            if A(j) == 1
                temp = M(j);
                samples = samples + 1;
                if temp == 1
                    occupied2 = occupied2 + 1;
                    n1 = n1 + 1;
                    % Max backoff set by user parameter
                    if n1 > q
                        n1 = q;
                    end
                    temp2 = j + n1 - 1;
                    if temp2 > length
                        temp2 = length;
                    end
                    A(j: temp2) = 0;
                    n2 = n2_def;
                elseif temp == 0
                    % Original algorithm----------------------------------
                    vacant2 = vacant2 + 1;
                    n1 = n1_def;
                    %-----------------------------------------------------
                    %  *************************************************
                    % Modified with restoration---------------------------
%                     vacant2 = vacant2 + 1;
%                     n2 = n2 + 1;
%                     if n2 > q
%                         n2 = q;
%                     end
%                     temp2 = j + n2 - 1;
%                     if temp2 > q
%                         temp2 = q;
%                     end
%                     A(j: temp2) = 1;      % Change from paper's algorithm
%                     n1 = n1_def;
                    %-----------------------------------------------------
                end 
            elseif A(j) == 0
                %n1 = n1_def;
            end                
        end

        % Calculate metrics     
        samplesTot(x, y) = samples;
        reductionTot(x, y) = samples / length;  
        vacanciesTot(x, y) = vacant;
        vacanciesTot2(x, y) = vacant2;
        occupTot(x, y) = occupied;
        occupTot2(x, y) = occupied2;
    end   
end

efficiency = 100.*(length - samplesTot)./(length);
vacanciesRatio = vacanciesTot2 ./ vacanciesTot;
vacanciesRatio(sweepsP, :) = 1;
opt = vacanciesRatio ./ reductionTot;