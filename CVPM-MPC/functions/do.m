function x=do(x,u,Rob,i,Ts, options)
% System
    tau=@(x)feedback_lin(x(1:2),x(3:4),u,Rob);
    dx=@(t,x)[x(3:4,1); forward_dyn(x(1:2), x(3:4,1), tau(x), Rob)];
    
    if i==1
        sol = ode23(dx, [0 Ts], x0, options);
    else
        sol = odextend(sol, dx, Ts*i);
    end
    x=sol.y(:,end);