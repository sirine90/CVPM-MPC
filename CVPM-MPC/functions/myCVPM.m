classdef myCVPM < handle
    % Class for a MPC-CVPM controller
    
    properties
        %% System
        % x_{k+1} = A*x_k + B*u_k
        A           % System Matrix
        B           % Input Matrix
        nx          % Number of state dimensions x_k
        nu          % Number of input dimensions u_k
        
        %% distrubance
        G           % Noise Matrix
        Sigmaw      % Covariance matrix of distrubance
        SigmaX      % Covariance matrix of state sequence
        %% Control Settings
        % MPC optimization
        % J = Sum_{k=0}^{N-1}(x_k'*Q*x_k +u_k'*R*u_k )+x_N'*S*x_N
        K           % LQ control law
        Nmpc        % MPC Horizon
        Ncvpm       % CVPM Horizon   
        %% Constraints
        constraint = struct(...
            'ineq',struct('E_',[],'e1',[],'e2',[],'F_',[],'f_',[]),...
            'eq',struct('E_',[],'e1',[],'e2',[],'F_',[],'f_',[]))  
        Pxt         % Polyhedron for terminal staes x_N
        Px0     % Feasible initial states
        
        Pu          % Polyhedron for inputs u_k for k in [0,N-1]
        Pw          % Polyhedron for noise w_k for k in [0,N-1]
        Pxp         % Probabilistic Constraint for the 4 states
        PXp         % Robust Probabilistic Constraint for complete State sequence of Ncvpm states i.e. (Xp^N \ominus G W^N)
        PT          % Target set for Case2 
        PU_ad       % Admissible Set for input sequence with Nmpc inptus
        PU_ad_cvpm  % Admissible Set for input sequence with Ncvpm inptus
        Pxp12       % Polyhedron for the first two states x_1 and x_2
        %% Precalculated data for MPC optimization
        A_mpc       % Lifted system matrix for Matrix eqution over whole Horizon
        B_mpc       % Lifted input matrix for Matrix eqution over whole Horizon
        H           % Hessian
        g           % Vector of Objective is f= g*x0
        Q_          % lifted weights
        %% Precalcualted data for CVPM optimizaton 
        A_cvpm      % Lifted system matrix for Ncvpm steps
        B_cvpm      % Lifted input matrix for Ncvpm steps
        G_cvpm      % Lifted disturbance matrix for Ncvpm steps
        %% Solver
        optionsMPC  % option for solver
        optionsCVPM % Option for Linear Program solver
        %% Logged Data
        data_x      % Log-matix for initial states
        data_X      % Log-matix for all state predictions
        data_U      % Log-matrix for input sequence
        data_case   % Log cases
        data_X_ref  % reference Trajectory 
        log_case    % Log-matrix for cases
        %% Plot
        fig
        u_max
        x_max
    end
    
    methods
        function obj = myCVPM(varargin)
            % Constructor for CVPM class
            % Inputs
            % 'A': (matix) system matrix
            % 'B': (matix) input matrix
            % 'G': (matix) distrubance matrix
            % 'Sigma': (matrix) covariance matrix of disturbance
            % 'Nmpc': (scalar) MPC horizon
            % 'Q': (matix) weight matrix for states
            % 'R': (matix) weight matrix for inputs
            % 'Pu': (Polyhedron) input set
            % 'Pw': (Polyhedron) disturbance set
            % 'Stability': (string, ['Lypunov','Riccati']) Mehtod for terminal cost
            % 'Ncvpm': (scalar) CVPM horizon
            
            %% Input Parser
            p = inputParser;
            p.addParameter('A',[]);
            p.addParameter('B',[]);
            p.addParameter('Nmpc',[]);
            p.addParameter('Q',[]);
            p.addParameter('R',[]);
            p.addParameter('Pu',[]);
            p.addParameter('Stability','Riccati');
            p.addParameter('Ncvpm', 1);
            p.addParameter('G', []);
            p.addParameter('Pw', []);
            p.addParameter('Sigma', []);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            %% input parser assignment
            obj.A=p.Results.A;
            obj.B=p.Results.B;
            obj.Nmpc=p.Results.Nmpc;
            obj.Pu=p.Results.Pu;
            % Dimensions
            obj.nx = size(obj.A,1);
            obj.nu = size(obj.B,2);
      
            %  Constraints
            if isempty(obj.Pu)
                obj.Pu = Polyhedron.fullSpace(obj.nu);
            end
            obj.Ncvpm =  p.Results.Ncvpm;
            obj.G =  p.Results.G;
            obj.Pw =  p.Results.Pw;
            obj.Sigmaw =  p.Results.Sigma;
            
            %% Init MPC
            % initFunction for MPC Optimizer
            % Terminal cost
            if strcmp( p.Results.Stability, 'Riccati')
                [Qt,~,obj.K] = dare(obj.A,obj.B,p.Results.Q,p.Results.R);
            elseif strcmp( p.Results.Stability, 'Lyapunov')
                Qt = dlyap(obj.A,p.Results.Q);
                obj.K = zeros(size(obj.B))';
            else
                error('Wrong value for Stability Option. Use Riccati or Lyapunov');
            end
            % _Precalc for minJ_
            % System Matrix for whole Horizon
            % X = A_*x0 + B_*U
            % X = [x1, x0, ... , xN]
            % U = [u0, u1, ... , u(N-1)]
            % A_
            obj.A_mpc =cell2mat(cellfun(@(x)obj.A^x,num2cell((1:obj.Nmpc)'),'UniformOutput',0));
            %  B_
            obj.B_mpc=tril(cell2mat(cellfun(@(x)obj.A^x,num2cell(toeplitz(0:obj.Nmpc-1)),'UniformOutput',0)))*kron(eye(obj.Nmpc),obj.B);
            % Hessian matrix and vecotor of Objective function
            %  J = 1/2 * U'*(B_'*Q_*B_ +R_)*U  + (B_'*Q_*A_*x0)' * U;
            %  J = 1/2 * U'*H*U  + x0' * g' * U;
            % H & g
            obj.Q_ = blkdiag(kron(eye(obj.Nmpc-1),p.Results.Q),Qt) ;
            obj.H = obj.B_mpc'* obj.Q_* obj.B_mpc +   kron(eye(obj.Nmpc),p.Results.R);
            obj.H =(obj.H+obj.H')/2; % ensure symmetric
            %obj.g = obj.B_mpc'*Q_*obj.A_mpc;

            %% Input sequence Set
            obj.PU_ad = cartesianProd(obj.Pu,obj.Nmpc);
            obj.PU_ad_cvpm = cartesianProd(obj.Pu,obj.Ncvpm);
            
            %% CVPM
            % initFunction for MPC Optimizer
            obj.A_cvpm =cell2mat(cellfun(@(x)obj.A^x,num2cell((1:obj.Ncvpm)'),'UniformOutput',0));
            obj.B_cvpm=tril(cell2mat(cellfun(@(x)obj.A^x,num2cell(toeplitz(0:obj.Ncvpm-1)),'UniformOutput',0)))*kron(eye(obj.Ncvpm),obj.B);
            obj.G_cvpm=tril(cell2mat(cellfun(@(x)obj.A^x,num2cell(toeplitz(0:obj.Ncvpm-1)),'UniformOutput',0)))*kron(eye(obj.Ncvpm),obj.G);
            % Covariance Matrix
            obj.calcSigma();
            % Initialize probabilistic Constraint to full space
            obj.setConstraint( Polyhedron.fullSpace(obj.nx) );
            %% Option for solver
            obj.optionsCVPM= optimoptions('quadprog','Display','off');
            obj.optionsMPC = optimoptions('quadprog','Display','off','MaxIter',1500,'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6);
            obj.initConstraints(varargin{:});
        end
              
        %% execute system prediction
        function Xn = F(obj,x0,U,n)
            % Calculate fist to nth system states depeding on input the
            % input sequence U
            if ~exist('n','var')
                n=1;
            elseif n*obj.nu > length(U)
                n = fix(length(U)/obj.nu);
            end
            Xn = obj.A_mpc(1:n*obj.nx,:) *x0 + obj.B_mpc(1:n*obj.nx,1: n*obj.nu) * U(1: n*obj.nu);
        end
   %% Reformulate constraints for optimization
        function initConstraints(obj,varargin)
            % Control Invariant Terminal Set
            obj.Pxt = ControlInvariantSet( obj.A, obj.B, obj.Pxp, obj.Pu );
            % State and Terminal Constraints depending on initial state:
            %   E_ *  U <= e1   - e2*x0
            %  Eeq_ * U  = e1eq - e2eq*x0
            obj.constraint.ineq.E_ = blkdiag(kron(eye(obj.Nmpc-1),obj.Pxp.A),obj.Pxt.A)*obj.B_mpc;
            obj.constraint.ineq.e1 = [kron(ones(obj.Nmpc-1,1),obj.Pxp.b);obj.Pxt.b];
            obj.constraint.ineq.e2 = blkdiag(kron(eye(obj.Nmpc-1),obj.Pxp.A),obj.Pxt.A)*obj.A_mpc;
            obj.constraint.eq.E_   = blkdiag(kron(eye(obj.Nmpc-1),obj.Pxp.Ae),obj.Pxt.Ae)*obj.B_mpc;
            obj.constraint.eq.e1   = [kron(ones(obj.Nmpc-1,1),obj.Pxp.be);obj.Pxt.be];
            obj.constraint.eq.e2   = blkdiag(kron(eye(obj.Nmpc-1),obj.Pxp.Ae),obj.Pxt.Ae)*obj.A_mpc;
            
            % Input Constraints:
            %   F_ * U <= f_
            %  Feq_ * U = feq_
            obj.constraint.ineq.F_ = kron(eye(obj.Nmpc),obj.Pu.A);
            obj.constraint.ineq.f_ = kron(ones(obj.Nmpc,1),obj.Pu.b);
            obj.constraint.eq.F_ = kron(eye(obj.Nmpc),obj.Pu.Ae);
            obj.constraint.eq.f_ = kron(ones(obj.Nmpc,1),obj.Pu.be);
            
            obj.Px0 = obj.calcFeasibleInitialStates();
        end
        
        %% calculate Feasible initial states
        function P = calcFeasibleInitialStates(obj)
            % Need first myMPC.init()
            E = obj.constraint.ineq.E_;
            e1 = obj.constraint.ineq.e1;
            e2 = obj.constraint.ineq.e2 ;
            F = obj.constraint.ineq.F_;
            f = obj.constraint.ineq.f_;
            E_eq = obj.constraint.eq.E_;
            e1_eq = obj.constraint.eq.e1;
            e2_eq = obj.constraint.eq.e2 ;
            F_eq = obj.constraint.eq.F_;
            f_eq = obj.constraint.eq.f_;
            
            % Create Polyhedron for all inputs U and init state x0
            P = Polyhedron('A',[e2 ,E;zeros(size(F,1),size(e2,2)), F],'b',[e1;f], ...
                'Ae',[e2_eq ,E_eq;zeros(size(F_eq,1),size(e2_eq,2)), F_eq],'be',[e1_eq;f_eq]);
            
            % Project constraints for all inputs U and init state x0 to the
            % intit state x0
            P.minHRep();
            P = P.projection(1:obj.nx);
        end        
%% Calculate admissible input set
        function P = calcAdmissibleInputs(obj,x0,varargin)
            
            if nargin == 3
                n = varargin{1};
            else
                n = 1;
            end
            E = obj.constraint.ineq.E_;
            e1 = obj.constraint.ineq.e1;
            e2 = obj.constraint.ineq.e2 ;
            F = obj.constraint.ineq.F_;
            f = obj.constraint.ineq.f_;
            E_eq = obj.constraint.eq.E_;
            e1_eq = obj.constraint.eq.e1;
            e2_eq = obj.constraint.eq.e2 ;
            F_eq = obj.constraint.eq.F_;
            f_eq = obj.constraint.eq.f_;
            
            % Create Polyhedron for all inputs U
            P = Polyhedron('A',[E;F],'b',[e1;f]-[e2;zeros(size(f,1),size(e2,2))]*x0,...
                'Ae',[E_eq;F_eq],'be',[e1_eq;f_eq]-[e2_eq;zeros(size(f_eq,1),size(e2_eq,2))]*x0);
            
            % Project constraints for all inputs U to the first n inputs+
            P.minHRep();
            if n<obj.Nmpc 
                P = P.projection(1:obj.nu*n);
            end
        end
        %% set probabilistic Constraint
        function setConstraint(obj,Pxp,transient)
            % determine the proabilistic constraint, such that the whole
            % state sequence do not violate it
            if ~exist('transient','var')
                transient='standard';
            end
            obj.Pxp=Pxp;
            obj.PXp = cartesianProd( obj.Pxp , obj.Ncvpm)- obj.G_cvpm * cartesianProd( obj.Pw , obj.Ncvpm);
            obj.PXp.minHRep;
            if obj.PXp.isEmptySet
                warning('PXp is calcFeasibleInitialStatesan Empty Set')
            end
            
            if strcmp(transient,'stable')
                S0 = obj.SigmaX(1:obj.nx,1:obj.nx)^-1;
                SN = obj.SigmaX((obj.Ncvpm-1)* obj.nx+1:end,(obj.Ncvpm-1)* obj.nx+1:end)^-1;
                % Initial States for Case 1:
                Px0_c1 = Polyhedron([obj.PXp.A*obj.A_cvpm,obj.PXp.A*obj.B_cvpm;zeros(size(obj.Pu.A,1)*obj.Ncvpm,obj.nx), kron(eye(obj.Ncvpm),obj.Pu.A)],[obj.PXp.b;kron(ones(obj.Ncvpm,1),obj.Pu.b)]);
                Px0_c1 = Px0_c1.projection(1:obj.nx);
                if Px0_c1.isEmptySet
                    error('Initial State for Case 1 does not exists: Stabiltiy is not possible')
                end
                Rw =  max(diag(obj.Pw.V*obj.G' *SN* obj.G*obj.Pw.V'));
                PXp2 =PolyMinusEllipse(Px0_c1,S0/Rw);
                obj.PT = cartesianProd( PXp2 , obj.Ncvpm);
                if obj.PT.isEmptySet
                    warning('The target set PT is an Empty Set')
                end
            else
                obj.PT =obj.PXp;
                %obj.PT =cartesianProd(obj.Pxp,obj.Ncvpm);
            end
        end
        
        %% Calculate Covariance matrix
        function calcSigma(obj,transient)
            % Determine the covariance matrix for the state sequence
            % regarding stability
            if ~exist('transient','var')
                transient='standard';
            end
            if strcmp(transient,'standard')
                obj.SigmaX = obj.G_cvpm*kron(eye(obj.Ncvpm),obj.Sigmaw)*obj.G_cvpm';
            elseif strcmp(transient,'stable')
                % Sigma for stability:
                obj.SigmaX =  kron(eye(obj.Ncvpm),obj.G'*obj.Sigmaw*obj.G);
                obj.SigmaX((obj.Ncvpm-1)*obj.nx+1:end,(obj.Ncvpm-1)*obj.nx+1:end) = dlyap(obj.A,(obj.G'*obj.Sigmaw*obj.G)^-1)^-1;
            else
                error('wrong option')
            end
        end
        
        %% execution of CVPM and MPC
        function u = do(obj,x0,x_ref)
            % performe MPC optimization and CVPM optimization
            % CVPM
            
            %option1: PU_ad_cvpm = cartesianProd(obj.Pu,obj.Ncvpm)
            [PU_cvpm, cvpm_case] = obj.CVPM(x0,obj.PU_ad_cvpm);
            
            
            % Optimization
            U = obj.MPC(x0,PU_cvpm,x_ref);
            u = U(1:obj.nu);
            % Data Logging
            obj.log_case    = [obj.log_case, cvpm_case];
            obj.data_x      = [obj.data_x,x0];
            obj.data_X      = [obj.data_X,obj.F(x0,U,obj.Nmpc)];
            obj.data_U      = [obj.data_U,U];
            obj.data_case   = [obj.data_case,cvpm_case];
        end
        
        %% MPC optimization function
        function U = MPC(obj,x0,varargin)
            % Need first myMPC.init()
            % U = argminJ(obj,x0)
            % U = argminJ(obj,x0,Polyhedron)
            
            % Constraints:
            
            Pmpc = obj.PU_ad;
            
            if nargin >= 3 % addtional constraints (e.g. CVPM set)
                Pmpc = Pmpc.intersect(varargin{1});
            end
            
            %Q_ = blkdiag(kron(eye(obj.Nmpc-1),p.Results.Q),Qt) ;
            %obj.g = obj.B_mpc'*Q_*obj.A_mpc;
            obj.g = obj.B_mpc'*obj.Q_*obj.A_mpc*x0 - obj.B_mpc'*obj.Q_*reshape(varargin{2},[],1);
            % Optimization:
            [U,~,exitflag,~] = quadprog(obj.H,(obj.g)',Pmpc.A,Pmpc.b,Pmpc.Ae,Pmpc.be,[],[],zeros(size(obj.H,1),1),obj.optionsMPC);
            
            
            % Error handling:
            if exitflag<=0
                warning('problem solving qp');
                disp(exitflag);
                U = zeros(obj.nu*obj.Nmpc,1);
            end
        end
        
        %%  CVPM optimization
        function [PU, cvpm_case] = CVPM(obj,x0,PU_ad)
            % Compute CVPM optimization
            % calc input set, that fulfills the probabilistic constraint
            
            PUp = invAffineMap2(obj.PXp,obj.B_cvpm,obj.A_cvpm*x0);
            PUp.minHRep;
            
            % Case Analysis
            PU = PUp.intersect(PU_ad);
            if ~PU.isEmptySet
                PU.minHRep;
                fprintf(' 1 ')
                cvpm_case = 1;
            else
                fprintf(' 2 ')
                cvpm_case = 2;
                

                %
                P= PU_ad * obj.PT;
                invSigma = obj.SigmaX^-1;
                H_cvpm = [obj.B_cvpm'*invSigma*obj.B_cvpm, -obj.B_cvpm'*invSigma;
                          -invSigma*obj.B_cvpm,                       invSigma];
                H_cvpm=(H_cvpm+H_cvpm')/2;
                % marginal weighting of the inputs for the analysis
                %%%H=(H+H')/2+1e-2* [eye(obj.nu*obj.Ncvpm),zeros(obj.nu*obj.Ncvpm,obj.nx*obj.Ncvpm);zeros(obj.nx*obj.Ncvpm,obj.nu*obj.Ncvpm),zeros(obj.nx*obj.Ncvpm)];
                f = x0'*obj.A_cvpm'* [invSigma*obj.B_cvpm,-invSigma ];
                %
                
                UX = quadprog(H_cvpm,f,P.A,P.b,P.Ae,P.be,[],[],[],obj.optionsCVPM);
                U_opt = UX(1:obj.nu*obj.Ncvpm);
                PU = Polyhedron(U_opt');
            end
            % Lift Polyhedron for all inputs
            PU = PU*Polyhedron.fullSpace(obj.nu*(obj.Nmpc-obj.Ncvpm));
            PU.minHRep;
        end
        
        %% Plot
        function subplot(obj,fig_nr)
            % Detailed plots of the simulation
            set(0,'DefaultLegendAutoUpdate','off')
            if  ~exist('fig_nr','var')
                fig_nr =1;
            end
            figure(fig_nr);
            clf;
            T = size(obj.data_x,2);
            % Max values
            obj.u_max = max( max(max(obj.Pu.V)), max(abs(obj.data_U(:))))*1.05;
            obj.x_max = max(abs(obj.data_x(:)))*1.05;
            %% States
            subplot(2,2,1)
            plot(0:T-1, obj.data_x(1:2,:), 0:T-1, obj.data_X_ref(1:2,1:T), 'LineWidth',2)
            if obj.nx==4
                legend('x_{1}','x_{2}', 'x_ref{1}', 'x_ref{2}')
            end
            hold on
            % case 1
            idx = find(obj.data_case==1)-1;
            plot(repmat(idx,2,1),repmat([-obj.x_max;obj.x_max],1,length(idx)),'color',[0.7,1,0.7])
            % case 2
            idx = find(obj.data_case==2)-1;
            plot(repmat(idx,2,1),repmat([-obj.x_max;obj.x_max],1,length(idx)),'color',[1,0.7,0.7])
            hold off
            axis([0 T -obj.x_max obj.x_max])
            xlabel('time')
            ylabel('States')
            title('States (Case 1: green; Case 2: red)')
            %% Inputs
            subplot(2,2,2)
            stairs(0:T-1, obj.data_U(1:obj.nu,:)','-o','LineWidth',2);
            if obj.nu==1
                legend('u')
            elseif obj.nu==2
                legend('u_1','u_2')
            elseif obj.nu==3
                legend('u_1','u_2','u_3')
            end
            hold on
            % case 1
            idx = find(obj.data_case==1)-1;
            plot(repmat(idx,2,1),repmat([-obj.u_max;obj.u_max],1,length(idx)),'color',[0.7,1,0.7])
            % case 2
            idx = find(obj.data_case==2)-1;
            plot(repmat(idx,2,1),repmat([-obj.u_max;obj.u_max],1,length(idx)),'color',[1,0.7,0.7])
            hold off
            axis([0 T -obj.u_max obj.u_max])
            xlabel('time')
            ylabel('Inputs')
            title('Inputs (Case 1: green; Case 2: red')
            %% Feasible Inital states
            subplot(2,2,3)
            if obj.nx == 4
                hold on
                plot(obj.data_x(1,:),obj.data_x(2,:),'r--o', obj.data_X_ref(1,1:T), obj.data_X_ref(2,1:T), 'r' )
                obj.Pxp12.plot('color','blue','alpha',0.3)
                %obj.Pxp.plot('color','blue','alpha',0.3)
                hold off
                axis equal
                axis([-obj.x_max obj.x_max -obj.x_max obj.x_max])
                xlabel('x_1/x_ref1')
                ylabel('x_2/x_ref2')
            end
            title('Feasible initial states(green) & CVPM(blue) & Trajectory(red)')
            %% Eigenvalues
            subplot(2,2,4)
            t = linspace(0,2*pi,90);
            unitcirc = cos(t)+1i*sin(t);
            plot(unitcirc);
            hold on
            Eig= eig(obj.A);
            plot(real(Eig),imag(Eig),'o','MarkerFaceColor','r')
            Eig= eig(obj.A-obj.B*obj.K);
            plot(real(Eig),imag(Eig),'o','MarkerFaceColor','y')
            hold off
            legend('Unit Circle','A','A-B*K')
            title('Eigenvalues')
            xlabel('Re')
            ylabel('Imag')
            axis equal
            grid on
        end
        
        function plot(obj,fig_nr)
            % animated plot of the simulation
            if  ~exist('fig_nr','var')
                fig_nr =1;
            end
            obj.fig=figure(fig_nr);
            clf;
            T = size(obj.data_x,2);
           
            % Max values
            obj.x_max = max(abs(obj.data_x(:)))*1.05;
            %% Feasible Inital states
            %hold on
            if obj.nx==4
                % Trajectory
                %plot(obj.data_X_ref(1,1:T), obj.data_X_ref(2,1:T), 'o');
                
                for i=1:T
                    hold on
                    if obj.log_case(1,i) == 1
                         plot(obj.data_x(1,i),obj.data_x(2,i),'r--o');
                    elseif obj.log_case(1,i) == 2
                         plot(obj.data_x(1,i),obj.data_x(2,i),'b--o');

                    end

                end
                % Single Point of Trajecotry: Slider
                h_slider = plot(obj.data_x(1,1),obj.data_x(2,1),'g--o','MarkerFaceColor','g','MarkerSize',5);
                %Prediction: Button
                X_pred =[obj.data_x(:,1) ,reshape(obj.data_X(:,1),[obj.nx, obj.Nmpc])];
                h_predict = plot(X_pred(1,:),X_pred(2,:),'g--o','MarkerFaceColor','g','MarkerSize',3.5);
                set(h_predict,'Visible','off');                
                plot(obj.data_X_ref(1,1:T), obj.data_X_ref(2,1:T), 'k','LineWidth', 1 )
                obj.Pxp12.plot('color','blue', 'alpha', 0.1)
              
                hold off
                axis equal
                axis([0 360 0 360])
                set(gcf,'color','w');
                xlabel('x_1/x_{ref,1}')
                ylabel('x_2/x_{ref,2}')
                h_title = title(['Feasible initial states(green) & Trajectory(red) & CVPM(blue), Act. Time Step: ',num2str(0) ]);
                grid on
            end
            %% Slider
            slider = uicontrol('Parent',obj.fig,'Style','slider','Position',[20,20,500,40],...
               'Value',1,'min',1, 'max',T,'SliderStep', [1/(T-1), 1/(T-1)]);
            slider.Callback = @(src,event)CallbackSlider(obj,src,h_slider,h_predict,h_title);
            %% Button Prediction
            button_predict = uicontrol('Parent',obj.fig,'Style','togglebutton','Position',[550,20,100,40],'String','Prediction');
            button_predict.Callback = @(src,event)CallbackButtonPrediction(obj,src.Value,h_predict);
        end
        %% Callbacks
        function CallbackSlider(obj,src,h_slider,h_predict,h_title)
            i = fix(src.Value);
            X_pred =[obj.data_x(:,i) ,reshape(obj.data_X(:,i),[obj.nx, obj.Nmpc])];
            
            if obj.nx==4
                set(h_slider,'XData',obj.data_x(1,i),'YData',obj.data_x(2,i));
                set(h_predict,'XData',X_pred(1,:),'YData',X_pred(2,:));
            elseif obj.nx==3
                set(h_slider,'XData',obj.data_x(1,i),'YData',obj.data_x(2,i),'ZData',obj.data_x(3,i));
                set(h_predict,'XData',X_pred(1,:),'YData',X_pred(2,:),'ZData',X_pred(3,:));
            end
            set(h_title,'String',['Feasible initial states(green) & Trajectory(red) & CVPM(blue), Act. Time Step: ',num2str(i-1) ])
            drawnow;
        end
        function CallbackButtonPrediction(~,Value,h_predict)
            if Value
                set(h_predict,'Visible','on');
            else
                set(h_predict,'Visible','off');
            end
        end
    end
end

