% [R,M,C,b] = fem_matrix(P,T,type,v,f) returns fem matrices according to
% 'type'. 'type' is a vector of characters which may contain 'r' (Rigidity
% Matrix), 'm' (Mass Matrix), 'c' (Convexity Matrix) and/or 'b' (right hand
% term). 'P', 'T' are standard points and triangles (tetrahedra) matrices.
% 'v' is the velocity field needed for 'M', whereas 'f' is the independent
% term, needed for 'b'. 

function [R,M,C,b] = fem_matrix(P,T,type,v,f)

N = size(P,1);
NT = size(T,1);
d = size(P,2);

Rigi = [zeros(1,d);eye(d)];
Masa = (ones(d+1)+eye(d+1))/12;
Conv = (ones(d+1)-eye(d+1))/2;
Conv = Conv(end:-1:1,:);

R = sparse(N,N);
M = sparse(N,N);
C = sparse(N,N);
b = sparse(N,1);
for r =1:NT
    nodos = T(r,:);
    puntos = P(nodos,:);
    H = [ones(1,d+1);puntos'];
    if ismember('r',type) || ismember('c',type);
        G=H\Rigi;
    end
    if ismember('r',type)
        Rloc = abs(det(H))*G*G'/prod(1:d);
        R(nodos,nodos) = R(nodos,nodos) + Rloc;
    end
    if ismember('m',type)
        Mloc = abs(det(H))*Masa/prod(1:d);
        M(nodos,nodos) = M(nodos,nodos) + Mloc;
    end
    if ismember('c',type) || ismember('f',type)
        medios = (puntos + puntos([2,3,1],:))/2;
    end
    if ismember('c',type)
        Cloc = abs(det(H))*Conv*v(medios)*G'/6;
        C(nodos,nodos) = C(nodos,nodos) + Cloc;
    end
    if ismember('f',type)
        bloc = abs(det(H))*Conv*f(medios)/6;
        b(nodos) = b(nodos) +bloc;
    end
end

