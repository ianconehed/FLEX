%% initialization
[v_it,v_kt] = deal(zeros(N,1)+v_rest); %membrane potentials
[R_yt,R_it,R_kt,sc_R] = deal(zeros(N,t_total/dt)); %rates (sc_R is spikes)
[s_yt,s_it,s_kt] = deal(zeros(N,1)); %activations
[g_Ey,g_Eyi,g_Ei,g_Ii,g_Ek] = deal(zeros(N,1)); %conductances
[t_ref_i,t_ref_k] = deal(zeros(N,t_total/dt + D)); %refractory periods
[T_ijp,T_ijd,del_W_ji,del_P_ik,del_W_in,del_T_ijp,del_T_ijd,H_p,H_d,H_p_scaled,H_d_scaled] = deal(zeros(N,N)); %synapse specific traces, weight update
[T_pt,T_dt,dopa,dopa_rec,dopa_rec_dynamic,dopa_ff] = deal(zeros(1,t_total/dt + dt)); %mean trace for population at time t
if l>1
    v_it(N-2*num_VTA+1:N) = v_it_init(N-2*num_VTA+1:N);
    v_kt(N-2*num_VTA+1:N) = v_kt_init(N-2*num_VTA+1:N);
    s_it(N-2*num_VTA+1:N) = s_it_init(N-2*num_VTA+1:N);
    s_kt(N-2*num_VTA+1:N) = s_kt_init(N-2*num_VTA+1:N);
    R_it(N-2*num_VTA+1:N,2+D-1) = R_it_init(N-2*num_VTA+1:N);
    R_kt(N-2*num_VTA+1:N,2+D-1) = R_kt_init(N-2*num_VTA+1:N);
    t_ref_i(N-2*num_VTA+1:N,2+D:2+D+10) = t_ref_i_init(N-2*num_VTA+1:N,:);
    t_ref_k(N-2*num_VTA+1:N,2+D:2+D+10) = t_ref_k_init(N-2*num_VTA+1:N,:);
end
%% input per trial
% if l < num_trials
%     t_miy = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,unit_stim,t_stim,t_reward,T_R,1,0);
% elseif l == num_trials
%     t_miy = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,unit_stim,t_stim,t_reward,T_R,0,0);
% end
if l == num_trials
    t_miy = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,5,t_stim,t_reward,T_R,1,0);
elseif mod(l,3) == 0
    t_miy = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,1,t_stim,t_reward,T_R,0,0);
elseif mod(l,3) == 1
    t_miy = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,3,t_stim,t_reward,T_R,0,0);
else
    t_miy = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,5,t_stim,t_reward,T_R,1,0);
end
%% learning per trial
%     eta_hebb = (l>num_trials/8).*eta_hebb_l;
%     eta_ff = (l>num_trials/8).*eta_ff_l;
%     eta_rec = (l>num_trials/8).*eta_rec_l;

eta_hebb = 0*(l>1).*eta_hebb_l;
eta_ff = 0*(l>1).*eta_ff_l;
eta_rec = 0*(l>1).*eta_rec_l;
%% time step loop
for t = 2+D:((t_total)/dt) %at each time step
    %% input conductance and FR
    s_yt = s_yt -(s_yt*dt/tau_si) + rho*(1-s_yt).*t_miy(:,t);
    R_yt(:,t) = R_yt(:,t-1) + (t_miy(:,t)/dt-R_yt(:,t-1))*(dt/tau_w);
    %% excitatory LIF
    exc_noise = [randn(N-num_VTA,1)/norm_noise;...
        randn(num_VTA,1)/norm_noise_VTA + 1.075*constant_VTA];

    del_vi = (exc_noise+g_L*(E_l-v_it) + (g_Ei + g_Ey)...
        .*(E_e - v_it) + g_Ii.*(E_i - v_it))*(dt/C_m);

    v_it = v_rest.*(t_ref_i(:,t)==1)...
        +(v_it + del_vi).*(v_it < v_th).*(t_ref_i(:,t)==0)...
        +v_hold.*(v_it>=v_th).*(t_ref_i(:,t)==0);

    sc_R(:,t) = (v_it == v_hold);
    t_ref_i(:,t:t+t_refractory) = t_ref_i(:,t:t+t_refractory) + (v_it == v_hold)*ones(1,t_refractory+1);
    s_it = s_it -(s_it*dt/tau_se) + rho*(1-s_it).*(v_it == v_hold);
    R_it(:,t) = R_it(:,t-1) + ((v_it == v_hold)/dt-R_it(:,t-1))*(dt/tau_w);
    %% inhibitory LIF
    inh_noise = [randn(N-num_VTA,1)/norm_noise;...
        randn(num_VTA,1)/norm_noise_VTA_inh + 1.1*constant_VTA];

    del_vk = (inh_noise +g_L*(E_l-v_kt) + (g_Ek + (iG/eG)*g_Eyi)...
        .*(E_e - v_kt)).*(dt./C_m_i);

    v_kt = v_rest.*(t_ref_k(:,t)==1)...
        +(v_kt + del_vk).*(v_kt < v_th_i).*(t_ref_k(:,t)==0)...
        +v_hold.*(v_kt>=v_th_i).*(t_ref_k(:,t)==0);

    t_ref_k(:,t:t+t_refractory) = t_ref_k(:,t:t+t_refractory) + (v_kt == v_hold)*ones(1,t_refractory+1);
    s_kt = s_kt -(s_kt*dt/tau_si) + rho*(1-s_kt).*(v_kt == v_hold);
    R_kt(:,t) = R_kt(:,t-1) + ((v_kt == v_hold)/dt-R_kt(:,t-1))*(dt/tau_w);

    %% Hebbian Activations
    hebb_rec = R_it(:,t-D)*R_it(:,t-dt)';
    hebb_ff = R_it(:,t-D)*R_yt(:,t-dt)';



    H_d(rec_identity) = eta_d*hebb_rec(rec_identity)/T_max_d;
    H_d(ff_identity) = eta_d1*hebb_ff(ff_identity)/T_max_d1;
    H_p(rec_identity) = eta_p*hebb_rec(rec_identity)/T_max_p;
    H_p(ff_identity) = eta_p1*hebb_ff(ff_identity)/T_max_p1;

    H_p_scaled(total_identity) = H_p(total_identity)/(1+100*dopa(t-1));
    H_d_scaled(total_identity) = H_d(total_identity)/(1+100*dopa(t-1));
    %% Traces
    del_T_ijp(rec_identity) = (-T_ijp(rec_identity) + H_p_scaled(rec_identity).*(T_max_p - T_ijp(rec_identity)))*(dt/tau_p);
    del_T_ijp(ff_identity) = (-T_ijp(ff_identity) + H_p(ff_identity).*(T_max_p1 - T_ijp(ff_identity)))*(dt/tau_p1);
    del_T_ijd(rec_identity) = (-T_ijd(rec_identity) + H_d_scaled(rec_identity).*(T_max_d - T_ijd(rec_identity)))*(dt/tau_d);
    del_T_ijd(ff_identity) = (-T_ijd(ff_identity) + H_d(ff_identity).*(T_max_d1 - T_ijd(ff_identity)))*(dt/tau_d1);

    T_ijp = max(T_ijp+del_T_ijp,0);
    T_ijd = max(T_ijd+del_T_ijd,0);
    %% Weight Updates
    del_W_ji(rec_identity) = del_W_ji(rec_identity)...
        + eta_rec*(T_ijp(rec_identity)-T_ijd(rec_identity))*(dopa(t-1)/num_VTA)*dt;

    if Ach(t-1)>dopa(t-1)
        Adp = dopa(t-1);
    else
        Adp = dopa(t-1);
    end
    del_W_in(ff_identity) = del_W_in(ff_identity)...
        + eta_ff(ff_identity).*(T_ijp(ff_identity)-T_ijd(ff_identity))...
        *(Adp/num_VTA)*dt;

    temp_hebb = R_it(:,t-D)*R_kt(:,t-dt)';
    %         temp_hebb = R_kt(:,t-D)*R_it(:,t-dt)';
    del_P_ik(hebb_identity) = del_P_ik(hebb_identity)...
        + eta_hebb(hebb_identity).*temp_hebb(hebb_identity)*(dopa(t-1))/num_VTA - lambda_P*P_ik(hebb_identity);

    %% Updating conductances
    g_Ey = W_in*s_yt; %input conductance
    g_Eyi = W_inin*s_yt; %input conductance to inhibitory cells
    g_Ei = W_ji*s_it; %recurrent excitatory conductance
    g_Ii = M_ki*s_kt; %I to E conductance
    g_Ek = P_ik*s_it; %E to I conductance

    %% dopamine calc
    temp_dopa = R_it(N-num_VTA+1:N,t-1)*1000;
    dopa(t) = dopa_func(mean(temp_dopa),thresh_size,d_0);
    if  l == 1
        Ach(t) = dopa(t);
    end
    %% save traces for plotting
    if select_rec == 1
        T_pt(t) = mean(T_ijp(1:npp,1:npp),'all')*100000; %recurrent
        T_dt(t) = mean(T_ijd(1:npp,1:npp),'all')*100000; %recurrent
    else
        T_pt(t) = mean(T_ijp(N-num_VTA+1:N,1:npp),'all')*5000000; %ff
        T_dt(t) = mean(T_ijd(N-num_VTA+1:N,1:npp),'all')*5000000; %ff
    end
end
%% update weights
W_ji = sparse(max(W_ji + del_W_ji,0));
W_in = sparse(max(W_in + del_W_in,0));
P_ik = sparse(max(P_ik + del_P_ik,0));
%% plotting variables
if plot_var == 1
    dopa_plot_pot = dopa;
    dopa_plot_dep = dopa;

    rec_vect(l) = mean(W_ji(1:npp,1:npp)*(npp^2),'all');
    m_vect(l) = mean(P_ik(N-num_VTA+1:N,npp+1:(2*npp))*(npp^2),'all');
    ff_vect(l) = mean(W_in(N-num_VTA+1:N,1:npp)*(npp^2),'all');
    if num_columns == 2
        rec_vect1(l) = mean(W_ji((2*npp+1):3*npp,(2*npp+1):3*npp)*(npp^2),'all');
        m_vect1(l) = mean(P_ik(N-num_VTA+1:N,(3*npp+1):(4*npp))*(npp^2),'all');
        ff_vect1(l) = mean(W_in(N-num_VTA+1:N,(npp+1):2*npp)*(npp^2),'all');
    end

    for t = 0:length(1:ROC_dt:t_total-2*ROC_dt)
        if l > 1
            hist_vect(:,t+1,l) = sum(sc_R(N-num_VTA+1:N,(t*ROC_dt+1):(t*ROC_dt + ROC_dt)),2);
        elseif l == 1
            hist_vect1(:,t+1) = sum(sc_R(N-num_VTA+1:N,(t*ROC_dt+1):(t*ROC_dt + ROC_dt)),2);
        end
    end
end
sc_R_it(:,:,l) = sc_R; %spiking for plotting
R_it_l(:,:,l) = R_it; %rates for plotting

sprintf('Trial %d complete',l)
%% trial plotting
if plot_var ==1
    if l == 1
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,num_VTA,R_it*1000,R_kt*1000,T_pt,T_dt,W_ji,'Before Learning',t_stim,dopa_plot_pot,dopa_plot_dep)
    elseif mod(l,1) == 0
        tit = sprintf('During Learning Trial %d',l);
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,num_VTA,R_it*1000,R_kt*1000,T_pt,T_dt,W_ji,tit,t_stim,dopa_plot_pot,dopa_plot_dep)
    elseif l == num_trials+1
        plot_func(T,delta,l,num_columns,0:dt:t_total,t_total,dt,npp,N,num_VTA,R_it*1000,R_kt*1000,T_pt,T_dt,W_ji,'After Learning (one stim)',t_stim,dopa_plot_pot,dopa_plot_dep)
    end
    drawnow
end
R_it_init = R_it(:,t_total);
R_kt_init = R_kt(:,t_total);
t_ref_i_init = t_ref_i(:,t_total:t_total+10);
t_ref_k_init = t_ref_k(:,t_total:t_total+10);
s_it_init = s_it;
s_kt_init = s_kt;
v_it_init = v_it;
v_kt_init = v_kt;
