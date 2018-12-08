% Main Program for the Geometry Block

% (C) 2018 Tiago Monteiro Cardoso
% Last update: 2018-12-08

clear
geometry_type = input ('Type of Geometry: ','s');
% m - monospheres; a - axons; b - axons 256
N =input ('Number of spheres/axons: ');
% up to 1000. In the case of spheres, 500 have been used; as for the axons, the respective submodule only outputs 256 axons.
PHI =input ('Fraction of Occupancy (packing): ');
% real between 0 and 1 - occupancy / degree of myelination
RES =input ('Voxel Resolution: ');
% so far, only used 256.
filename =input ('File name: ','s');
% including the extension ('.mat')

if geometry_type == 'm'
    Diam = 2*RES*(((3*PHI)/(4*pi*N))^(1/3));

    connected_matrix = burning_algo( uint32(monospheres( N,PHI,RES )),RES );

    save(filename, 'connected_matrix', 'Diam') 
end

if geometry_type == 'a'
    
    Diam = 10;
    
    connected_matrix = uint32(axons(N,PHI,RES));
    
    save(filename, 'connected_matrix', 'Diam')
end

if geometry_type == 'b'
    
    Diam = 10;
    
    connected_matrix = uint32(axons_256(N,PHI,RES));
    
    save(filename, 'connected_matrix', 'Diam')
end



