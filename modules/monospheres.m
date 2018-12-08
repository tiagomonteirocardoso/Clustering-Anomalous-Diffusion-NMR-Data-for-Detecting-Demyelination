% Submodule for the creation of simulation spaces consisting of spheres of equal diameter.

% (C) 2018 Tiago Monteiro Cardoso
% Last update: 2018-12-08


function bin_matrix = monospheres( N,PHI,RES )
% Creates monospheres

%Creates a voxel filled with spheres

% Calculated Radius 
r = RES*(((3*PHI)/(4*pi*N))^(1/3));

% Initialization of the centers' coordinates
X0 = (RES-2*r)*rand(N,3)+r*ones(N,3);

% Bounds
lb = zeros(N,3)+r*ones(N,3);
ub = RES*ones(N,3)-r*ones(N,3);

% Unused constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Function Potential
function f = mono_sphere_potential( C )
% calculates the potential in the case of non-penetrating mono spheres

S = size(C);
num_spheres = S(1);
f = 0;

for i=1:num_spheres
    for j=1:num_spheres
        if i < j
% The best up to now for higher concentrations
            f = f + (10000)*((1 - ((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)/((2*r)^2))^4)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r));
% 100 spheres 0.5 / 20 spheres 0.6 / 100 spheres 0.6
%             f = f + (10000)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r))^4)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r));
% 20 spheres 0.5
%             f = f + (10000)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r))^2)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r));
% Original
%             f = f + (1/2)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r))^2)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2 + (C(i,3)-C(j,3))^2)^(1/2))/(2*r));
        end
    end
end

end

% Minimization Procedure (= final coordinates of the spheres' centers)
X = fmincon(@mono_sphere_potential,X0,A,b,Aeq,beq,lb,ub);

% Binarization of the volumes
x=0.5:1:(RES-0.5);
y=x;
z=x;
bin_matrix=zeros(RES,RES,RES);
o=1;

for n=1:N
    for k=1:RES
        for l=1:RES
            for m=1:RES
                if (x(k)-X(n,1))^2+(y(l)-X(n,2))^2+(z(m)-X(n,3))^2 <= r^2
                    bin_matrix(k,l,m)=o;
                end
            end
        end
    end
    o = o + 1;
end


end

