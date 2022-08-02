function tau=feedback_lin(q,qd,u,Rob)
   tau=Rob.itorque(q', u')+ qd'*Rob.coriolis(q', qd') +Rob.gravload(q');
   tau=tau';