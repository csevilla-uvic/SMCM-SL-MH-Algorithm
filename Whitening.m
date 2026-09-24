%CS This will compute the Whitening Matrix
 
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
N = 2304; %CS Number of Sample Paths Per Distribution
alpha = 1; %CS U_rng keep proportion for proposal 
N_keep = round(alpha * N);

%CS My Vectors
Sample = zeros(L, dim);
phi = zeros(L, dim); %CS Our log-transform variable
Likeli = zeros(L, 1);

%CS Initial Sample For each Timescale
%CS I'm cheating here a bit
Sample(1, :) = [1.5, 2, 2, 3, 8.344, 8.541, 10];

%CS Log-Transform Save First Initial Step
phi(1, : ) = log(Sample(1, :));
S_x = zeros(N, 2*dim);    %CS Our Count Summary Statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Define my Prior and Proposal Distributions in Phi Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Prior Parameters in Phi Space, Note transformation to Log-Norm
%CS Prior Chosen to be mildly Informative
mu_prior = [2,2,2,2,2,2,2]; %Define in Phi Space (Normal)
sigma_prior = [0.5^2,0.5^2,0.5^2,0.5^2,0.5^2,0.5^2,0.5^2]; %Define in Phi Space (Normal), variances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Scan TXT documents for Data Set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Read Matrix for Predictive Variables
%CS Only U needed Now
U = readmatrix("U_11520_5min (2).txt",'Delimiter',{','});

%CS Coursening Step
%U = U(1:9:end, :);
U = U(T_start:T_end, :);

CAPEs = U(:,1);
CAPEls = U(:,2);
Drys = U(:,3);
Cns = U(:,4);
Wn = U(:, 5);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CS Here I'll code the MCMC Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%CS For Correlated Sample Paths We define our Seeds USE FOR every
u_rng = randi(1e9, N, 1);

%CS Compute Synthetic Likelihood for Initialization
Timescales = Sample(1, :);
Initial = [96,3,0,1];
for i=1:N
    %CS This Sample Path uses the Injection to compute Sx
    stream = RandStream('mrg32k3a', 'Seed', u_rng(i));    %CS Set seed for our substream
    [S_x(i, :)] = Sample_Path(Timescales, Initial, DT0, n, Niter, dim, CAPEs, CAPEls, Drys, Cns, Wn, stream);
end

%CS Compute the Multivariate Normal Likelihood. Assume Meanfield for now
old_mean = mean(S_x, 1);
old_sigma = cov(S_x);

%CS Find the Whitening Matrix
eigvals = eig(old_sigma);
epsilon = 5e-4 * median(eigvals);
old_sigma_reg = old_sigma + epsilon * eye(2*dim);

% Step D: whitening matrix
[U,D] = eig(old_sigma_reg);
W =  diag(1 ./ sqrt(diag(D))) * U';

% Step E: whiten observed summaries
%S_y_reg = (W * (S_y - old_mean)')';
%S_x_reg = (W * (S_x - old_mean)')';
%TEST = cov(S_x_reg);

fid1 = fopen('Whiten.txt','w');
fprintf(fid1, '%f', W);


