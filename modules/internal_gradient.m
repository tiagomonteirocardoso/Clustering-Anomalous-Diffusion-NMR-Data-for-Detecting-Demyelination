% Submodule for the computation of the field due to internal (susceptibility induced) gradients.
% It is necessary to manually 'activate' either line 12 or 15 (variable "logicalmat"), if the simulation being performed is 
% in a space containing respectively spheres or axons.

% (C) 2018 Tiago Monteiro Cardoso
% Last update: 2018-12-08

function B = internal_gradient(x,y,z,file_geometry,deltachi,B0)

load(file_geometry, 'connected_matrix', 'Diam')
% For the monospheres
logicalmat = (connected_matrix<=1000);

% For the axons
% logicalmat = (connected_matrix==1);
clear connected_matrix

B = 0;
% direction of B0 (must be kept this way. If changed, the expression for B
% should change)
direction = [1,0,0];
sizelogicalmat = size(logicalmat);

parfor i=1:9
    for j=1:9
        for k=1:9
%             rx = x-floor(x)+4.5-i;
%             ry = y-floor(y)+4.5-j;
%             rz = z-floor(z)+4.5-k;
%             r = sqrt(rx.^2 + ry.^2 + rz.^2);
%             cos_ang = (direction(1)*rx + direction(2)*ry + direction(3)*rz)./(sqrt(rx.^2 + ry.^2 + rz.^2));
            B = B + ((deltachi*(3*(((direction(1)*(x-floor(x)+4.5-i))./(sqrt((x-floor(x)+4.5-i).^2 + (y-floor(y)+4.5-j).^2 + (z-floor(z)+4.5-k).^2))).^2)-1).*B0)./(4*pi*((sqrt((x-floor(x)+4.5-i).^2 + (y-floor(y)+4.5-j).^2 + (z-floor(z)+4.5-k).^2)).^3))).*logicalmat(sub2ind(sizelogicalmat,floor(x)+i-4,floor(y)+j-4,floor(z)+k-4));
        end
    end
end
end
