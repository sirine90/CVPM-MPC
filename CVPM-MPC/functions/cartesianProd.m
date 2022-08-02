function [ OutPolyhedron ] = cartesianProd( InPolyhedron , n )
% Make a cartesian product: 
% P_in * P_in * P_in * .... = \prod_i^n P_in = P_in^n
%
% Inputs:
% InPolyhedron: (Polyhedron) input polyhedron, which is stacked n times
% n: (scalar) number of stack operations
%
% Outputs:
% OutPolyhedron: (Polyhedron) stacked polyhedron, Dimension is n times dimension of InPolyhedron

if isempty(InPolyhedron.He)
    OutPolyhedron = Polyhedron(kron(eye(n),InPolyhedron.A),kron(ones(n,1),InPolyhedron.b));
else
    OutPolyhedron = Polyhedron('A',kron(eye(n),InPolyhedron.A),'b',kron(ones(n,1),InPolyhedron.b),...
        'Ae',kron(eye(n),InPolyhedron.Ae),'be',kron(ones(n,1),InPolyhedron.be));
end
OutPolyhedron.minHRep;
end

