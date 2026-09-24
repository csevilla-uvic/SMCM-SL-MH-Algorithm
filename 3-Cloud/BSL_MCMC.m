%CS This Will be the main function for My Synthetic Liklielihood 
%CS Start by Reading a Data Set
 
clear all

%CS Set RNG for Whole Program
masterSeed = 12345;
stream = RandStream('mrg32k3a', 'Seed', masterSeed); %CS Creates re-usable seed Fixes the Starting point of any random number generator 
RandStream.setGlobalStream(stream); %CS Makes this seed global to code

%Note that we still have our RNG split into substreams
%Each of these is independent and reproducible

DT0 = 0.0833;
T_start = 1;   %CS Truncate our time series
T_end = 1440;
Niter = T_end+1 - T_start;   %CS Length of Time Series We want to use
n = 10; % system dimension
dim =7;
L = 25000; %CS Number of MCMC Samples
N = 192; %CS Number of Sample Pathsfor Sx
M = 192; %CS Number of Sample Paths for Sy
eps = 1e-6; %CS Regularization Constant 
alpha = 1; %CS U_rng keep proportion for proposal 
N_keep = round(alpha * N);

%CS My Vectors
Sample = zeros(L, dim);
phi = zeros(L, dim); %CS Our log-transform variable
Likeli = zeros(L, 1);

%CS Initial Sample For each Timescale
Sample(1, :) = [1,1,1,1,1,1,1];

%CS Log-Transform Save First Initial Step
phi(1, : ) = log(Sample(1, :));
S_x = zeros(N, 2*dim);    %CS Our Count Summary Statistics
S_y = zeros(M, 2*dim);    %CS Our Count Summary Statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Define my Prior and Proposal Distributions in Phi Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Prior Parameters in Phi Space, Note transformation to Log-Norm
%CS Prior Chosen to be mildly Informative
mu_prior = [0,0,0,0,0,0,0]; %Define in Phi Space (Normal)
sigma_prior = [10^2,10^2,10^2,10^2,10^2,10^2,10^2]; %Define in Phi Space (Normal), variances

sigma_prop = [0.5^2,0.5^2,0.5^2,0.5^2,0.5^2,0.5^2,0.5^2]; %CS Proposal Distribution Variances in terms of tau);
sigma_prior = diag(sigma_prior);

%CS Define the Whitening Matrix 
S_0 = [0.0795268494405864,0.0565360363618827,0.00997058256172841,0.0474684727044753,0.0664595992476852,0.0174765504436728,0.0461434823495370,0.336643034408198,0.269994357306172,0.101445492707848,0.225798081478528,0.259398240722777,0.133438273351449,0.215443235198209];   %CS Summary statistic mean used to find Whitening MAtrix

White = [56.7631362353083,-856.177432658182,-1175.44999435823,1093.91207303415,-58.0379812143334,1121.68753002542,-216.608960430957,1.52949349317981,3.48038459476295,43.3838507846397,-15.4714668590265,1.13639884167825,-63.0607313093159,6.64777101382701;
109.592027767395,-75.0274326092151,-520.359091221173,-701.403718430557,-124.026123754710,235.553840045537,879.712481457768,-2.20514917209598,4.53303043658931,58.4328763925226,24.4160685524719,8.63488429036015,-40.7619879049281,-68.0062466949652;
297.647857016732,162.893897453194,-703.542931332392,56.7080541032042,-343.507256605988,-743.311072872183,-240.814711237212,-2.87130944753700,-5.00309800629618,86.5304712414808,-8.71063314945336,21.5155077675106,140.105080969157,24.7984589515412;
303.542157578402,343.665420370539,42.0196012502738,-183.135933177750,-378.983797517699,385.523229566288,-289.563743835681,-4.03439081574645,-16.9010948987473,-4.95185508077267,13.8632435945860,39.1386969872932,-153.300492824422,54.6883181780205;
273.712224826975,-259.199914148307,264.896857557742,190.304948547212,-371.199710030821,-87.0035610059602,228.589957894951,-9.07285134437207,20.3065196516538,-109.433807923213,-23.9916851503085,52.9362012623384,63.8228951106496,-62.5403613921265;
-27.0301957457851,-289.015644070436,14.6336466882616,-226.479755735896,-42.0851889706080,-17.5610251884594,-182.078393911109,1.38462852118427,44.8368408040196,-21.3127792893234,64.8959028073147,25.6329632159897,37.9437569238440,148.240260301468;
-245.788298436248,28.5824052880168,-8.54211938783305,17.3254878958752,-176.647354738184,0.607342867563196,13.2478797885485,37.1682627259767,-7.77916529301877,23.2175189401772,-7.36847345591118,140.150880277038,-2.01523610026855,-10.7336708835280;
0.619614704909873,15.0707977173470,0.213629591095541,13.4333760713616,1.14554559873016,2.59652595507239,26.1284128164067,1.22992841028140,-21.5054791025069,-9.61210370692351,-109.112427630097,-0.0198633929142758,13.3051490735842,132.199033832257;
6.14386844587074,-17.1723627551882,4.25339590760023,-4.91000129876880,3.82336352329006,-11.8655869728021,-6.96107605505245,-11.9937866300852,75.3588475068931,50.0391252397498,-60.8513880320246,13.7954249022580,-65.4622187212754,-23.6285254602473;
24.0767562579641,2.75473716179684,2.30614953109414,2.93685182869723,27.5871173342251,2.58019367455756,2.75365557875881,-59.8482735308691,-15.0998278694722,19.9442564927745,14.2258567120805,89.2549600639022,13.7124984283624,8.36239671709707;
0.414476101433249,-0.881895099104553,6.77401967233450,-3.79859348106401,-5.87212579809100,9.68071173520845,-3.80852943047305,-1.22877668426945,-1.91682293284834,43.4645659637198,-11.6536077854991,-13.0405381446731,42.2112080571895,-9.92176392763343;
-0.245567869689575,-5.72833924977402,5.31573382568810,6.47997372117776,-5.43940328227553,-6.80141391704065,6.33423832295659,-1.81084218245631,-25.2324000259845,33.3403146910077,20.0029776852293,-11.6941878154575,-28.7179771465477,16.6505594277523;
7.49822333431039,4.12841755849536,1.33000482436319,4.41143121102488,5.89769389986049,0.942611366532110,4.32676225075064,28.9982986642738,15.0198977268429,8.01042953587944,12.5196479476683,12.2254441553061,4.23675005081258,10.7837297385874;
-5.25119529900437,6.37231565857672,-0.0902620229682786,5.09480381351821,-4.89577321624118,1.12674244132331,4.96327554648741,-17.9773791887306,21.2557707913099,-0.343070378706747,13.7750781032906,-10.0505447689899,4.55721248444602,12.3293176926532];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Scan TXT documents for Data Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Read Matrix for Predictive Variables
%CS Only U needed Now
U = readmatrix("U_11520_5min (2).txt",'Delimiter',{','});

%CS Coursening Step
%U = U(1:9:end, :);
U = U(T_start:T_end, :);
Niter = size(U,1);

CAPEs = U(:,1);
CAPEls = U(:,2);
Drys = U(:,3);
Cns = U(:,4);
Wn = U(:, 5);

%CS For Correlated Sample Paths We define our Seeds USE FOR every
u_rng = randi(1e9, N, 1);

%Perform Preprocessing By taking Sample Path Average
%I'm Cheating With the Initial And Timescales Here
Timescales = [1.5, 2, 16, 3, 8.344, 8.541, 10];
Initial = [96,3,0,1];
for i=1:M
    %CS This Sample Path uses the Injection to compute Sx
    stream = RandStream('mrg32k3a', 'Seed', u_rng(i));    %CS Set seed for our substream
    [S_y(i, :)] = Sample_Path(Timescales, Initial, DT0, n, Niter, dim, CAPEs, CAPEls, Drys, Cns, Wn, stream);
end
S_y = mean(S_y , 1);

%CS Use whitening Matrix
%S_y = (White * (S_y - S_0)')';


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Here I'll code the MCMC Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CS Open File For continuous Printing
fid2 = fopen('Sample.txt','w');
fprintf(fid2, '%f,', Sample(1, :));
fprintf(fid2, '\n');

%CS For Correlated Sample Paths We define our Seeds USE FOR every
u_rng = randi(1e9, N, 1);

%CS Compute Synthetic Likelihood for Initialization
Timescales = Sample(1, :);
for i=1:N
    %CS This Sample Path uses the Injection to compute Sx
    stream = RandStream('mrg32k3a', 'Seed', u_rng(i));    %CS Set seed for our substream
    [S_x(i, :)] = Sample_Path(Timescales, Initial, DT0, n, Niter, dim, CAPEs, CAPEls, Drys, Cns, Wn, stream);
end

%CS Compute the Multivariate Normal Likelihood. Assume Meanfield for now
%S_x = (White * (S_x - S_0)')';
%CS Let's Try some slight regularization for T=8640
old_mean = mean(S_x, 1);
old_sigma = cov(S_x);
old_sigma = old_sigma + eps*eye(size(old_sigma));
%[old_sigma, lambda] = cov1Para(S_x);   %CS Variance Reduction Using Ledoit & Wolf

pos_def = isequal(old_sigma,old_sigma') && all( eig(old_sigma)>0 ) && (cond(old_sigma)<1e+10);
if pos_def
    oldlikeli = logmvnpdf(S_y,old_mean,old_sigma);
else %CS If can't be resolved, automatically reject
    disp('Removed')
    disp(Timescales)
    oldlikeli = -1E99;
end
Likeli(1) = oldlikeli;
oldprior = logmvnpdf(phi(1, :), mu_prior, sigma_prior);

fid1 = fopen('Likeli.txt','w');
fprintf(fid1, '%f', Likeli(1, :));
fprintf(fid1, '\n');

for l=1:L 
    Sample(l+1,:) = Sample(l,:); 
    phi(l+1,:) = phi(l,:);
    for d=1:dim
        %CS Sample but in phi Space
        phi(l+1, d) = phi(l, d) + mvnrnd(0, sigma_prop(d));
    
        %CS New u_rng Proposals
        keep_idx = randperm(N, N_keep); 
        replace_idx = setdiff(1:N, keep_idx);
        new_seeds = randi(1e9, length(replace_idx), 1);
        u_proposed = u_rng;
        u_proposed(replace_idx) = new_seeds;
    
        %CS Find Our Sample Paths
        Timescales = exp(phi(l+1, :));
        for i=1:N
            stream = RandStream('mrg32k3a', 'Seed', u_proposed(i));    %CS Set seed for our substream
            [S_x(i, :)] = Sample_Path(Timescales, Initial, DT0, n, Niter, dim, CAPEs, CAPEls, Drys, Cns, Wn, stream);
        end
        %CS In Least Approach, Sometimes S_y can be zero, which is NO GOOD.
        %S_x(find(any(S_x == 0, 2))', :) = [];
        %S_x = (White * (S_x - S_0)')';
        new_mean = mean(S_x, 1);
        new_sigma = cov(S_x);
        new_sigma = new_sigma + eps*eye(size(new_sigma));
        %CS Need to Check Covariance is Positive Definite
        pos_def = isequal(new_sigma,new_sigma') && all( eig(new_sigma)>0 ) && (cond(new_sigma)<1e+10);
        if pos_def
            newlikeli = logmvnpdf(S_y,new_mean,new_sigma);
        else %CS If can't be resolved, automatically reject
            disp('Removed')
            disp(Timescales)
            newlikeli = -1E99;
        end
    
        newprior = logmvnpdf(phi(l+1, :), mu_prior, sigma_prior);
    
        accept = newlikeli + newprior  - oldlikeli - oldprior;
        accept = min(1, exp(accept));
        
        %CS Sometimes the log gaussian is positive 
        if unifrnd(0,1) <= accept
            oldlikeli = newlikeli;
            oldprior = newprior;
            u_rng = u_proposed;     %CS update rng Streams
            Sample(l+1, :) = Timescales;
        else 
            phi(l+1,d) = phi(l,d);
        end
    
    end
    Likeli(l+1) = oldlikeli;
    fprintf(fid2, '%f,', Sample(l+1, :));
    fprintf(fid2, '\n');
    fprintf(fid1, '%f', Likeli(l+1));
    fprintf(fid1, '\n');
end


