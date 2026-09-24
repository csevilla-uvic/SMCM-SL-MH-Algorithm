%CS This Will be the main function for My Synthetic Liklielihood 
%CS Start by Reading a Data Set
 
clear all

%CS Set RNG for Whole Program
masterSeed = 12345;
stream = RandStream('mrg32k3a', 'Seed', masterSeed); %CS Creates re-usable seed Fixes the Starting point of any random number generator 
RandStream.setGlobalStream(stream); %CS Makes this seed global to code

%Note that we still have our RNG split into substreams
%Each of these is independent and reproducible

DT0 = 0.2;
T_start = 1;   %CS Truncate our time series
T_end = 1440;
Niter = T_end+1 - T_start;   %CS Length of Time Series We want to use
n = 20; % system dimension
dim =7;
L = 25000; %CS Number of MCMC Samples
N = 192; %CS Number of Sample Paths Per Distribution
del_t = 0.02; %CS Proposal Distribution Sigma
alpha = 1; %CS U_rng keep proportion for proposal 
N_keep = round(alpha * N);

%CS Here I can include the gaps in data I want to ignore
skip = [];%[438:703, 854:1227, 1561:1682, 5737:5890]; %Int to skip
T_length = Niter - length(skip); %True length of data after skipiping

%CS My Vectors
Sample = zeros(L, dim);
phi = zeros(L, dim); %CS Our log-transform variable
Likeli = zeros(L, 1);

%CS Initial Sample For each Timescale
Sample(1, :) = [10,10,10,10,10,10,10];

%CS Log-Transform Save First Initial Step
phi(1, : ) = log(Sample(1, :));
S_x = zeros(N, 2*dim);    %CS Our Count Summary Statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Define my Prior and Proposal Distributions in Phi Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Prior Parameters in Phi Space, Note transformation to Log-Norm
mu_prior = [0,0,0,0,0,0,0]; %Define in Phi Space (Normal)
sigma_prior = [1,1,1,1,1,1,1]; %Define in Phi Space (Normal)

sigma_prop = diag(del_t*ones(1, dim));
sigma_prior = diag(sigma_prior);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Scan TXT documents for Data Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Read Matrix for Exogenous Variables
%CAPEs = readmatrix("CAPE.txt");
CAPEls = readmatrix("CAPE_Orig.txt");
%Drys = readmatrix("DRY.txt");
Cns = readmatrix("CIN_Orig.txt");
%Wn = readmatrix("Wn.txt");

%CAPEs = CAPEs(T_start:T_end);
CAPEls = CAPEls(T_start:T_end);
%Drys = Drys(T_start:T_end);
Cns = Cns(T_start:T_end);
%Wn = Wn(T_start:T_end);

%CS Read CLoud Data Matrix
Nc = readmatrix("Nc.txt", 'Delimiter', ',');
Nd = readmatrix("Nd.txt", 'Delimiter', ',');
Ns = readmatrix("Ns.txt", 'Delimiter', ',');

N_fc = Nc(T_start:T_end);
N_fd = Nd(T_start:T_end);
N_fs = Ns(T_start:T_end);
N_fcs = n*n - (N_fc + N_fd + N_fs);

%CS Create a vector of Number of Transitions at Each time step
Delta_Ncs = N_fcs(2:Niter) - N_fcs(1:Niter-1);
Delta_Nc = N_fc(2:Niter) - N_fc(1:Niter-1);
Delta_Nd = N_fd(2:Niter) - N_fd(1:Niter-1);
Delta_Ns = N_fs(2:Niter) - N_fs(1:Niter-1);

W = [Delta_Ncs, Delta_Nc, Delta_Nd, Delta_Ns];

%CS Run our Injection for the Data, only once is necessary for 3 clouds
Ncount = zeros(T_length, dim);

Initial = [N_fcs(1), N_fc(1), N_fd(1), N_fs(1)];

for I = 1:Niter-1
    if any(skip == I)  
    else
        K=sum(abs(W(I, :)))/2;
        if K > 4   %CS Seems like the fastest 
            [Ncount(I, :)]= echelon(W(I, :), dim, K, stream);
        else
            [Ncount(I, :)]= injection(W(I, :), dim, stream);
        end
    end
end
N_avg = sum(Ncount)/T_length;
N_var = std(Ncount);
S_y = [N_avg, N_var];


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Here I'll code the MCMC Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%CS For Correlated Sample Paths We define our Seeds USE FOR every
u_rng = randi(1e9, N, 1);

%CS Compute Synthetic Likelihood for Initialization
Timescales = Sample(1, :);
for i=1:N
    %CS This Sample Path uses the Injection to compute Sx
    stream = RandStream('mrg32k3a', 'Seed', u_rng(i));    %CS Set seed for our substream
    [S_x(i, :)] = Sample_Path(Timescales, Initial, DT0, n, Niter, T_length, skip, dim, CAPEls, Cns, stream);
end

%CS Compute the Multivariate Normal Likelihood. Assume Meanfield for now
old_mean = mean(S_x, 1);
old_sigma = cov(S_x);

pos_def = isequal(old_sigma,old_sigma') && all( eig(old_sigma)>0 ) && (cond(old_sigma)<1e+10);
if pos_def
    oldlikeli = logmvnpdf(S_y,old_mean,old_sigma);
else %CS If can't be resolved, automatically reject
    disp('Removed')
    disp(Timescales)
    oldlikeli = -1E99;
end

fid1 = fopen('Likeli.txt','w');
%fprintf(fid1, '%f,', norm(old_mean - S_y, 2));
fprintf(fid1, '%f,', oldlikeli);

%[old_sigma, lambda] = cov1Para(S_x);   %CS Variance Reduction Using Ledoit & Wolf

%pos_def = isequal(old_sigma,old_sigma') && all( eig(old_sigma)>0 ) && (cond(old_sigma)<1e+10);
%if pos_def
%    oldlikeli = logmvnpdf(S_y,old_mean,old_sigma);
%else %CS If can't be resolved, automatically reject
%    disp('Removed')
%    disp(Timescales)
%    oldlikeli = -1E99;
%end

%fprintf(fid1, '%f,', oldlikeli);

%old_sigma = diag(diag(old_sigma));

%pos_def = isequal(old_sigma,old_sigma') && all( eig(old_sigma)>0 ) && (cond(old_sigma)<1e+10);
%if pos_def
%    oldlikeli = logmvnpdf(S_y,old_mean,old_sigma);
%else %CS If can't be resolved, automatically reject
%    disp('Removed')
%    disp(Timescales)
%    oldlikeli = -1E99;
%end

%fprintf(fid1, '%f', oldlikeli);
%fprintf(fid1, '\n');

for l=1:L 
    u_rng = randi(1e9, N, 1);
    for i=1:N
        stream = RandStream('mrg32k3a', 'Seed', u_rng(i));    %CS Set seed for our substream
        [S_x(i, :)] = Sample_Path(Timescales, Initial, DT0, n, Niter, T_length, skip, dim, CAPEls, Cns, stream);
    end
    %CS In Least Approach, Sometimes S_y can be zero, which is NO GOOD.
    %S_x(find(any(S_x == 0, 2))', :) = [];

    new_mean = mean(S_x, 1);
%    fprintf(fid1, '%f,', norm(S_y - new_mean, 2));

    new_sigma = cov(S_x);

    %CS Need to Check Covariance is Positive Definite
    pos_def = isequal(new_sigma,new_sigma') && all( eig(new_sigma)>0 ) && (cond(new_sigma)<1e+10);
    if pos_def
        newlikeli = logmvnpdf(S_y,new_mean,new_sigma);
    else %CS If can't be resolved, automatically reject
        disp('Removed')
        disp(Timescales)
        newlikeli = -1E99;
    end

    fprintf(fid1, '%f,', newlikeli);

%    [new_sigma, lambda] = cov1Para(S_x);

        %CS Need to Check Covariance is Positive Definite
%    pos_def = isequal(new_sigma,new_sigma') && all( eig(new_sigma)>0 ) && (cond(new_sigma)<1e+10);
%    if pos_def
%        newlikeli = logmvnpdf(S_y,new_mean,new_sigma);
%    else %CS If can't be resolved, automatically reject
%        disp('Removed')
%        disp(Timescales)
%        newlikeli = -1E99;
%    end

%    fprintf(fid1, '%f', newlikeli);

%    new_sigma = diag(diag(new_sigma));

%    pos_def = isequal(new_sigma,new_sigma') && all( eig(new_sigma)>0 ) && (cond(new_sigma)<1e+10);
%    if pos_def
%        newlikeli = logmvnpdf(S_y,new_mean,new_sigma);
%    else %CS If can't be resolved, automatically reject
%        disp('Removed')
%        disp(Timescales)
%        newlikeli = -1E99;
%    end

%    fprintf(fid1, '%f', newlikeli);
%    fprintf(fid1, '\n');

end


