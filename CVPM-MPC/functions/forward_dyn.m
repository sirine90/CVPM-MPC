function ddq=forward_dyn(q,dq,tau,Rob)
     ddq=Rob.accel(q',dq',tau');