%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       DELACHEV ADDS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Function where I can change Timescales for SL
function     [fsc,fc,fd,fs, Delta_N] = SMC_only(Timescales,dim,fsc,fc,fd,fs,C,Cl,D,Cn,Wn,DT,n,stream)
%pass additional area fraction and extra largescale variables (or declare them as global)
%add additional area fractions to solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       DELACHEV ADDS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Number of sites = n^2
nsite = n^2;
%CS Save the Number of transitions that occur
Delta_N = zeros(1, dim);


% Conditional rates
[r01,r02,r03,r12,r13,r23,r34,r10,r20,r30,r40]=  equilibriumdistributionSL(Timescales,C,Cl,D,Cn,Wn); %IG additional rates, additional LS variables
% Define cloud state
cloud = [fsc,fc,fd,fs]; %IG additional area fracs
%min(cloud)

% define the seven possible transitions
diffv(1,1:4)  = [1,0,0,0];  %01
diffv(2,1:4)  = [0,1,0,0];  %02
diffv(3,1:4)  = [0,0,1,0];  %03
diffv(4,1:4)  = [-1,1,0,0]; %12
diffv(5,1:4)  = [-1,0,1,0]; %13
diffv(6,1:4)  = [0,-1,1,0]; %23
diffv(7,1:4)  = [0,0,-1,1]; %34
diffv(8,1:4)  = [-1,0,0,0]; %10
diffv(9,1:4)  = [0,-1,0,0]; %20
diffv(10,1:4) = [0,0,-1,0]; %30
diffv(11,1:4) = [0,0,0,-1]; %40
%IG add additional transition defs

% We will use Gilespie's algorithm and iterate until we have exceeded the
% time step, DT

time = 0;
count = 0;

num_trans = 0;

while time < DT
    count = count+1;
    % Calculate the number of clearsky elements
    clearsky = nsite - sum(cloud);

    % Absolute rates 
    R01 = r01*clearsky;
    R02 = r02*clearsky;
    R03 = r03*clearsky;
    R12 = r12*cloud(1);
    R13 = r13*cloud(1);
    R23 = r23*cloud(2);
    R34 = r34*cloud(3);
    R10 = r10*cloud(1);
    R20 = r20*cloud(2);
    R30 = r30*cloud(3);
    R40 = r40*cloud(4);

    %IG add additional rates
    
    Rmat = [R01,R02,R03,R12,R13,R23,R34,R10,R20,R30,R40]; %IG add additonal rates
    
    Rsum = cumsum(Rmat); %cumulative sum
    
    time;
    if(Rsum(dim) < 0) %IG Change to last index. (11)
        'RSUM ERROR'
            Rmat
             Rsum
             count
             time
             cloud
    
        pause
    elseif(Rsum(dim) == 0)  %IG Change to last index. (11)
        t = 2*DT;
    elseif(Rsum(dim) > 0) %IG Change to last index. (11)
        t = -1*log(rand(stream,1,1))/Rsum(dim); %IG Change to last index. (11)
    end    
    
    % Calculate the time until the next transition

    time = time+t;
    
%%reject transition if time exceeds DT 
    if(time>DT)
        break
    end
    
    % Which transition occurs?
    test = rand(stream, 1,1);
    Rsum = Rsum/Rsum(dim); %IG Change to last index. (11)
    Rsum = Rsum - test;

    tmp =  find(sign(Rsum)==1);
    
    ind = min(tmp);

    diff = sum(diffv(ind,:),1);
    Delta_N(ind) = Delta_N(ind) + 1;
    %IG not sure what's happening in the above 3 lines. prob don't need to change anything but try to make sense.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       DELACHEV ADDS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_trans = num_trans + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       DELACHEV ADDS        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


    % New cloud array
    cloud = cloud + diff;
    
    
end
fsc = cloud(1);
fc = cloud(2);
fd = cloud(3);
fs = cloud(4);

%IG add additonal area fraction.

%num_trans

