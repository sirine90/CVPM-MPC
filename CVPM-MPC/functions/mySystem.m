classdef mySystem < handle
    % Class for the propatation of the system with different distrubance
    % options
    
    properties
        Rob     % Robot Manipulator
        G       % Uncertainty Matrix
        Px      % polyhedral state set
        Pu      % polyhedral input set
        Pw      % polyhedral distrubance set
        Sigma   % Covariance matrix for distrubance
        Ts      % Sampling Time
        option  % kind of distrubance use
        sol
    end
    
    methods
        %% Constructor
        function obj = mySystem(Rob, G, Px, Pu, Pw, Sigma, Ts, option)
            % Constructor
            % Inputs:
            % Px: (Polyhedron) state set
            % Pu: (Polyhedron) input set
            % G: (matrix) disturbance matrix
            % Pw: (Polyhedron) distrubance set
            % Sigma: (matrix) covariance matrix of distrubance
            % options: (string)
            %           'none': no distrubance
            %           'uni': use Polyhedron.randomPoint for distrubance
            %           'gauss': disturbance is gaussian distributed 
             

            if nargin == 5
                obj.Rob=Rob;
                obj.Px=Px;
                obj.Pu=Pu;
                obj.Pw=Pw;
                obj.G=G;
                obj.Ts=Ts;
                obj.Sigma=[];
                obj.option='uni';
            elseif nargin == 6
                obj.Rob=Rob;
                obj.Px=Px;
                obj.Pu=Pu;
                obj.Pw=Pw;
                obj.G=G;
                obj.Ts=Ts;
                obj.Sigma=Sigma;
                obj.option='none';
            elseif nargin == 7
                obj.Rob=Rob;
                obj.Px=Px;
                obj.Pu=Pu;
                obj.Pw=Pw;
                obj.G=G;
                obj.Ts=Ts;
                obj.Sigma=Sigma;
                obj.option='gauss';
            elseif nargin == 8
                obj.Rob=Rob;
                obj.Px=Px;
                obj.Pu=Pu;
                obj.Pw=Pw;
                obj.G=G;
                obj.Ts=Ts;
                obj.Sigma=Sigma;
                obj.option=option;
            else 
                error('wrong input')
            end
        end
        
        %% execute one step of the sysem
        function x=do(obj,x,u,i)
            % Execution of one step of the system, starting at state x0 and
            % using the input u. The distribution selected in the constructor 
            % is used
            tau=@(x)feedback_lin(x(1:2),x(3:4),u,obj.Rob);
            dx=@(t,x)[x(3:4,1); forward_dyn(x(1:2), x(3:4,1), tau(x), obj.Rob)];
            
            options = odeset('RelTol',1e-5);
            if i==1
                obj.sol = ode23(dx, [0 obj.Ts], x, options);
            else
                obj.sol = odextend(obj.sol, dx, obj.Ts*i);
            end
            x=obj.sol.y(:,end);
       
            if strcmp(obj.option,'none')
                x=x;
            elseif strcmp(obj.option,'uni')
                %% TODO: randomPoint is not uniform
                x= x + obj.G*obj.Pw.randomPoint;
            elseif strcmp(obj.option,'gauss')
                w = mvnrnd(zeros(obj.Pw.Dim,1),obj.Sigma,1)';
                while ~obj.Pw.contains(w)
                    w = mvnrnd(zeros(obj.Pw.Dim,1),obj.Sigma,1)';
                end
                x= x + obj.G*w;
            else
                x = x +obj.G*obj.Pw.randomPoint;
            end
            % Check bounds
            if ~obj.Pu.contains(u)
                warning('System violate input constraints')
            elseif  ~obj.Px.contains(x)
                warning('System violate state constraints')
            end
               
        end
    end
end

