function y = t_mik1(n,num_VTA,dt,T,t_total,p_r,npp,unit_stim,t_stim,t_reward,T_R,alpha,beta)
t_mik11 = zeros(n,t_total/dt);
ndd = npp;
for k = unit_stim
    j = (((k+1)/2)-1)*ndd;
    t_mik11(j+1:j+ndd,(t_stim((k+1)/2)+dt)/dt:(t_stim((k+1)/2)+T)/dt) = (rand(ndd,T/dt) < p_r);
%     t_mik1(j+npp+1:j+2*npp,1:t_stim/dt) = (rand(npp,t_stim/dt) < .5*p_r);
%     t_mik11(n-num_VTA+1:n,(t_stim(k)+dt)/dt:(t_stim(k)+T)/dt) = (rand(num_VTA,T/dt) < p_r);
end
if T_R ~=0
    t_mik11(n-num_VTA+1:n,(t_reward+dt)/dt:(t_reward+T_R)/dt) = (rand(num_VTA,T_R/dt) < alpha*p_r);
    if beta ==1
        t_mik11(n-num_VTA+1:n,(t_stim(2)+dt)/dt:(t_stim(2)+T_R)/dt) = (rand(num_VTA,T_R/dt) < alpha*p_r);
    end
end
y = t_mik11;
end