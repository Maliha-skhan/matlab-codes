clear all;close all; clc;

G=[1,0;0,0.7] % rank is 2

%%%%%% THE GOOD %%%%%%

%Assume (unrealistically!) that we have measured noise-free data:
d_pure=[0.5;0.001]

%Find the solution ?_pure to the problem, using simple matrix algebra
m_pure = (G)^(-1)*d_pure

%Assume now that there is noise on the data:
n=[0.008;0.011]

% the norm of n
n_norm=norm(n)

% What is the signal-to-noise ratio (ratio of norms)
s_N_ratio= norm(d_pure)/(n_norm)

% Find a solution ? to the inverse problem for the noise contaminated data
m_noise = (G)^(-1)*(d_pure+n)
s_Nm_ratio= norm(m_pure)/norm(m_noise-m_pure) % signal to noise ratio of the noisy m vector


%%%%%% THE BAD %%%%%%
Gb=[1,0.002;0,0] % rank is 1

%How many solutions to (1) exist in this case (? = ??) ?
% infinetly many solutions
% anything can be put in 

%mb= [0.5; alpha]% alpha is any real number %(Gb)^(-1)*d_pure


%%%%%% THE UGLY %%%%%%
Gu=[1,0;0.002,10^(-24)] % rank is 2
m_u = (Gu)^(-1)*d_pure
m_unoise = (Gu)^(-1)*(d_pure+n) % the 1st comp is 0.508 but MATLAB is rounding to zero
s_N_u= norm(m_u)/norm(m_unoise-m_u)  %signal to noise ratio % the solution is completely dead 

