function x_ref = sys_trajectory(x0,v,T,steps,N)
T_sim=T;
qref_final=[50 50];
dqref_final=[0 0];
[q_ref, dq_ref, ddq_ref] = jtraj(x0(1:2)', qref_final, linspace(0,T_sim+1,steps+N));
x_ref= [q_ref' ; dq_ref'];
end

