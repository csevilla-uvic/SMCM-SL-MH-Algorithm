%This script to USE the Synthetic Likelihood
%CS NOTE DIFFERENT FROM OTHER USES

function [S_x] = Sample_Path(Timescales, Initial, DT0, n, Niter, dim, CAPE, CAPEL, DRY, CIN, Wn, stream)
  
CAPE0 = 200;
N2=0.0001; %s^-2, Vaisala frequency, squared,
ZT=16000;% height of the troposphere
c=sqrt(N2)*ZT/pi;% ~50 m/sec; gravity wave speed; velocity time scale
CAPE0 = CAPE0/c^2;
Ntrans = zeros(Niter, dim);

N_fcs = zeros(1, Niter+1);
N_fsc = zeros(1, Niter+1);
N_fc = zeros(1, Niter+1);
N_fd = zeros(1, Niter+1);
N_fs = zeros(1, Niter+1);

N_fcs(1) = Initial(1);
N_fsc(1) = Initial(2);
N_fc(1) = Initial(3);
N_fd(1) = Initial(4);
N_fs(1) = Initial(5);

%CS Save Transitions at each time Step
%CS All randomness must come from rand / randn (global stream)
for I=1:Niter
    [N_fsc(I+1),N_fc(I+1),N_fd(I+1),N_fs(I+1), Delta_N]=SMC_only(Timescales,dim,N_fsc(I),N_fc(I),N_fd(I),N_fs(I),CAPE(I)/CAPE0,CAPEL(I)/CAPE0, DRY(I), CIN(I), Wn(I),DT0,n, stream);
    N_fcs(I+1) = n*n - (N_fsc(I+1) + N_fc(I+1) + N_fd(I+1) + N_fs(I+1));
    Ntrans(I, :) = Delta_N;
end
%CS Create a vector of Number of Transitions at Each time step
Delta_Ncs = N_fcs(2:Niter+1) - N_fcs(1:Niter);
Delta_Nsc = N_fsc(2:Niter+1) - N_fsc(1:Niter);
Delta_Nc = N_fc(2:Niter+1) - N_fc(1:Niter);
Delta_Nd = N_fd(2:Niter+1) - N_fd(1:Niter);
Delta_Ns = N_fs(2:Niter+1) - N_fs(1:Niter);

W = [Delta_Ncs', Delta_Nsc', Delta_Nc', Delta_Nd', Delta_Ns'];


%CS Adding our New Minimum Cost Flow Algorithm
Ncount = flow(W')';

N_avg = sum(Ncount)/Niter;
N_var = std(Ncount);
S_x = [N_avg, N_var];

