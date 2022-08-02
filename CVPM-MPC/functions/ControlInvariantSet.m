function [Xf, isMax] = ControlInvariantSet( A, B, X, U , varargin)
% Computes the maximum Control Invariant Set starting from X
%
% Inputs
% A: System matrix
% B: Input matrix
% X: Polyhedron for x
% U: Polyhedron for u
% option: Procciding if algroithm dosen not converge:
%         approx: Use a control invriant Set as approximation for the
%         maximum control invariant set
%         origin: Use origin as  terminal constraint
%         else:   Use outcome of the abort calculation
%
% Outpus
% Xf: control invariant set
% isMax: (Boolean) true if method is converaged


validationScalar = @(x) validateattributes(x,{'numeric'},{'scalar'});
validationLogic = @(x) validateattributes(x,{'logical'},{'scalar'});
p = inputParser;
addParameter(p,'option','approx');
addParameter(p,'maxItr',20,validationScalar);
addParameter(p,'maxItrApprox',5,validationScalar);
addParameter(p,'print',true,validationLogic);
p.parse(varargin{:});

maxItr       = p.Results.maxItr;
maxItrApprox = p.Results.maxItrApprox;
option       = p.Results.option;
print        = p.Results.print;
converged    = false;

for i = 1:maxItr
    Xf = PreSet(A,B, X, U);
    Xf = Xf.intersect(X).minHRep();
    if Xf==X
        converged = true;
        isMax=true;
        break
    else
        X = Xf;
    end
end
if print
    fprintf('Iterations for maximum control invariant set: %d\n', i);
end

if ~converged
    isMax = false;
    if strcmp(option,'approx')
        if print
            fprintf('Computation of maximum control invariant set finished without convergence. \nTherefore use a approximated control invariant set.\n');
        end
        Xf = Polyhedron('Ae',eye(size(A)),'be',zeros(size(A,1),1)); % only origin
        for i = 1:maxItrApprox
            Xf2 = PreSet(A,B, Xf, U);
            Xf2=Xf2.intersect(X);
            if Xf2==Xf
                converged = true;
                break
            else
                Xf = Xf2;
            end
        end
        Xf = Xf.intersect(X).minHRep();
        if print
            fprintf('Iterations for control invariant set: %d\n', i);
        end
    elseif strcmp(option,'origin')
        if print
            fprintf('Computation of maximum control invariant set finished without convergence. \nTherefore use origin.\n');
        end
        Xf = Polyhedron('Ae',eye(size(A)),'be',zeros(size(A,1),1)); % only origin
    else
        if print
            fprintf('Computation of maximum control invariant set finished without convergence.');
        end
    end
end

end

