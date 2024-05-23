%% structure parameters
num_columns = 1;
npp = 100;
ppc = 2*npp;
num_VTA = 100;
N = ppc*num_columns + num_VTA;
pop = num_columns*2;
n = N;
unit_num = 1:2:(2*num_columns);
% num_trials = 600;
num_trials = 80;

%% time parameters
dt = 1; %time step
T = 100/dt; %time stimulated
T_R = 150/dt; %time of reward
delta = 600/dt; %time in between stimulation
t_total = 1801; %time of trial
D = 10/dt; %intrinsic delay
tau_w = 40; %estimate time window
tau_dopa = 40;
tau_max = 50;

%% input parameters
eG = .04; %input gain excitatory
iG = .04; %input gain inhibitory
VTA_g = .005;
W_in = input_weights(num_columns,npp,N,num_VTA,eG,VTA_g);
W_inin = (iG/eG)*W_in;
W_inin(N-num_VTA:N,:) = 0;
unit_stim = unit_num;
t_stim = [201 701 1101 2201 4001];
t_reward = 1100;
p_r = .03*dt; %poisson rate 40 Hz

%% weight matrix initialization
l5_rec = .00014;%.00016; DO NOT SET TO ZERO OR RECURRENT IDENTITY WILL BE WRONG
l23_rec = .00000;
l5_l23 = .0005; %.0002
l23_l5 = .0000002;%.4; %DO NOT SET TO ZERO OR FF IDENTITY WILL BE WRONG
i_l5_rec = .0001;
i_l23_rec = .000;
i_l5_l23 = .02;%.07;%.05
i_l23_l5 = 0;
i_l23_l23 = .0;
i_l5_l5 = .0;
l5_VTA = 0.00001;
p_scale = -.1;
p_l5_rec = .0003;
p_l23_rec = .001;
p_l5_l23 = 0;
p_l23_l5 = 0;
l_23_l_VTA = .00000000000000000000001;
l_VTA_l_VTA = .0015;

W_ji = Sparse_L_ij(num_columns,npp,N,num_VTA,l5_rec,l23_rec,l5_l23,l23_l5,0,0,0,0,.000);
M_ki = Sparse_L_ij(num_columns,npp,N,num_VTA,i_l5_rec,i_l23_rec,i_l5_l23,i_l23_l5,i_l23_l23,i_l5_l5,0,l_VTA_l_VTA,.0000);
P_ik = Sparse_L_ij(num_columns,npp,N,num_VTA,p_l5_rec,p_l23_rec,p_l5_l23,p_l23_l5,0,0,l_23_l_VTA,0,.0001);
rec_identity = W_ji.*L_ij_no_rand(num_columns,npp,N,num_VTA,1,0,0,0,0,0,0)>0;
ff_identity = W_in.*L_ij_no_rand(num_columns,npp,N,num_VTA,0,0,0,0,0,1,0)>0;
hebb_identity = P_ik.*L_ij_no_rand(num_columns,npp,N,num_VTA,0,0,0,0,0,0,1)>0;
total_identity = sparse((rec_identity + ff_identity)>0);

%% membrane dynamics parameters
rho = 1/7; %percentage change of synaptic activation with input spikes
tau_se = 80; %time constant for excitatory synaptic activation
tau_si = 20; %time constant for inhibitory synaptic activation
tau_si_VTA = 5;
tau_se_VTA = 20;
norm_noise = 1e5;
norm_noise_VTA = 20;
constant_VTA = .035;
norm_noise_VTA_inh = 10;
C_m = .2; %membrane capacitance
C_m_i = C_m*ones(N,1);
C_m_i(N-num_VTA+1:N) = .1;
g_L = .01; %leak conductance
E_i = -70;
E_l = -60; %leak reversal potential
E_e = -5; %excitatory reversal potential
v_th = -55; %threshold potential
v_th_i = -50; 
v_rest = -60; %resting potential
v_hold = -61; %return potential
t_refractory = 2/dt;

%% Learning parameters
tau_p = 1800; %LTP time constant
tau_d = 800; %LTD time constant
T_max_p = .003; %maximum eligibility trace for LTP
T_max_d = .0033; %maximum eligibility trace for LTD
eta_p = 300; %
eta_d = 135; %

tau_p1 = 2000;% tau_p1 = 200; %LTP time constant
tau_d1 = 800;% tau_d1 = 800; %LTD time constant
T_max_p1 = .0015;% T_max_p1 = .003; %maximum eligibility trace for LTP
T_max_d1 = .004;% T_max_d1 = .0034; %maximum eligibility trace for LTD
eta_p1 = 650;% eta_p1 = 25*3000; %
eta_d1 = 40;% eta_d1 = 15*3000; %
trace_refractory = zeros(N,N);


eta_rec_l = .00015; %learning rate
eta_hebb_l = [zeros(N)];
eta_hebb_l(201:300,101:200) = abs(.1*rand(npp,npp));
eta_ff_l = [abs(1.25*rand(N,N))];
eta_p_rec = 1;
K = 1;
select_rec = 1;

num_thresh_rec = 20;
num_thresh_ff = 0;


%% additional parameters

ROC_dt = 50;
[ff_vect, ff_vect1] = deal(zeros(1,num_trials+1));
[rec_vect, rec_vect1] = deal(zeros(1,num_trials+1));
[m_vect, m_vect1] = deal(zeros(1,num_trials+1));
hist_vect = zeros(length(N-num_VTA+1:N),(t_total-1)/ROC_dt,num_trials);
hist_vect1 = zeros(length(N-num_VTA+1:N),(t_total-1)/ROC_dt);
auc_plot = zeros(length(N-num_VTA+1:N),(t_total-1)/ROC_dt);
lambda = 0;
temp_t = zeros(num_VTA,t_total/dt);
rew_cons = 2;

%% activation function parameters
thresh_size = 4;
d_0 = 5;

%% balance parameters

h_0 = 0.005;
tau_adap = 4000;
adap_act = 4;
h_max = 2;