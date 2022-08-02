function [PolyOut] = PolyMinusEllipse(PolyIn,Emat)
% Pontriagin difference of a Polyhedron
% PolyIn  = {x| PolyIn.A * x <= PolyIn.b}
% with a eliptic set E = {x|x'*Emat*x <= 1}
% yields
% PolyOut = {x| PolyIn.A * x <= PolyIn.b - f}
% with f_i = max_{x in E} PolyIn.A(i,:)*x
%
% Inputs
% PolyIn: (Polyhedron) Polyhedron, which is the subtrahend of the Pontriagin difference 
% Emat: (matrix) this matrix defines the ellipse E = {x|x'*Emat*x <= 1}, which is the minuend of the pontriagin difference
%
% Outputs
% PoylOut: (Polyhedron) pontiagin difference, represented as polyhedron

It = size(PolyIn.A,1); % Number of inequalitys
n = size(PolyIn.A,2);  % Dimenstion of space
f = zeros(It,1);       % Initialize of offset vector f
    function [c,ceq]=ellipsoid(x,Emat) % Function for ellipsoids
        c = x'*Emat*x-1;
        ceq = 0;
    end
ellipsoid2 = @(x)ellipsoid(x,Emat); % anonymous funciton of  ellipsoid with Matrix Emat

options = optimoptions('fmincon','Display','off');
%sequence of linear programs with ellipitc constraints
% TODO: find better solver than fmincon
for i = 1:It
    f(i) =PolyIn.A(i,:)* fmincon(@(x)-PolyIn.A(i,:)*x,zeros(n,1),[],[],[],[],[],[],ellipsoid2,options);
end

% Output Polyhedron
PolyOut = Polyhedron(PolyIn.A,PolyIn.b-f);

end

