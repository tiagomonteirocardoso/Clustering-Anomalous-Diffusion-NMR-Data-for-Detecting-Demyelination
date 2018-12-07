function random_walk_function_axons_II
% same as Random_walk but it's a function, unique file and for the axons

filename = 'interaxon_post256_PHI_80_DEM_60_XIII.mat';

% num_spins = 3300;
num_spins = 1000;

% step_meter = (10^(-6))*18/50;
step_meter = (10^(-6))*(0.1);

D = (2.3)*(10)^(-9);
TE = (84.4)*(10)^(-3);
% TE = (44.4)*(10)^(-3);

load(filename, 'connected_matrix')

RES = length(connected_matrix);
dt = ((step_meter)^2)/(6*D);
numsteps = floor(TE/dt);

% 17.6 (PHI = 80(entire voxel)/88%(selected part of voxel)); 16 (PHI =
% 70/80%); 12.65 (80% in the selected part of post. body); 8.8 (80%
% Splenium 256); 6.7 (80% Posterior Body 256); 7.8 (70% Splenium 256)

convert_coord_meter = (6.7)*(10^6);
% convert_coord_meter = (8.8)*(10^6);
% convert_coord_meter = (7.8)*(10^6);
% convert_coord_meter = 16*(10^6);
% convert_coord_meter = (12.65)*(10^6);

step = step_meter*convert_coord_meter;


X = zeros(numsteps + 1,num_spins,3);
X0 = 96*rand(num_spins,3)+80;
% X0 = 176*rand(num_spins,3)+40;
% X0 = 216*rand(num_spins,3)+20;

for i1 = 1:num_spins
    continue1 = 1;
    while continue1
        a1 = floor(X0(i1,1))+1;
        b1 = floor(X0(i1,2))+1;
        c1 = floor(X0(i1,3))+1;
        
%         if a1==RES+1
%             a1=RES;
%         end
%         if b1==RES+1
%             b1=RES;
%         end
%         if c1==RES+1
%             c1=RES;
%         end
        
        if connected_matrix(a1,b1,c1) == 0
            continue1 = 0;
        else
%             X0(i1,1) = 96*rand+80;
            X0(i1,1) = 96*rand+80;
            X0(i1,2) = 96*rand+80;
            X0(i1,3) = 96*rand+80;
        end
    end
end

X(1,:,:) = X0;

for i2 = 2:numsteps+1
    for j2 = 1:num_spins
        continue2 = 1;
        i3 = 1;
        
        while continue2
            
            crossover = 0;
            
            X(i2,j2,1) = X(i2-1,j2,1) + (step/sqrt(3))*randn;
            X(i2,j2,2) = X(i2-1,j2,2) + (step/sqrt(3))*randn;
            X(i2,j2,3) = X(i2-1,j2,3) + (step/sqrt(3))*randn;
                
            
            if X(i2,j2,1) < 20
                crossover = 1;
                X(i2,j2,1) = X(i2,j2,1) + (fix((20-X(i2,j2,1))/216)+1)*216;
            elseif X(i2,j2,1) >= 236
                crossover = 1;
                X(i2,j2,1) = X(i2,j2,1) - fix((X(i2,j2,1)-20)/216)*216;
            end
            
            if X(i2,j2,2) < 20
                crossover = 1;
                X(i2,j2,2) = X(i2,j2,2) + (fix((20-X(i2,j2,2))/216)+1)*216;
            elseif X(i2,j2,2) >= 236
                crossover = 1;
                X(i2,j2,2) = X(i2,j2,2) - fix((X(i2,j2,2)-20)/216)*216;
            end
            
            if X(i2,j2,3) < 20
                crossover = 1;
                X(i2,j2,3) = X(i2,j2,3) + (fix((20-X(i2,j2,3))/216)+1)*216;
            elseif X(i2,j2,3) >= 236
                crossover = 1;
                X(i2,j2,3) = X(i2,j2,3) - fix((X(i2,j2,3)-20)/216)*216;
            end
            
            a2 = floor(X(i2,j2,1))+1;
            b2 = floor(X(i2,j2,2))+1;
            c2 = floor(X(i2,j2,3))+1;
            
            a_3_4 = floor((X(i2,j2,1)+(X(i2,j2,1)+X(i2-1,j2,1))/2)/2)+1;
            b_3_4 = floor((X(i2,j2,2)+(X(i2,j2,2)+X(i2-1,j2,2))/2)/2)+1;
            c_3_4 = floor((X(i2,j2,3)+(X(i2,j2,3)+X(i2-1,j2,3))/2)/2)+1;
            
            a_medio = floor(abs((X(i2,j2,1)+X(i2-1,j2,1))/2))+1;
            b_medio = floor(abs((X(i2,j2,2)+X(i2-1,j2,2))/2))+1;
            c_medio = floor(abs((X(i2,j2,3)+X(i2-1,j2,3))/2))+1;
            
            a_1_4 = floor((X(i2-1,j2,1)+(X(i2,j2,1)+X(i2-1,j2,1))/2)/2)+1;
            b_1_4 = floor((X(i2-1,j2,2)+(X(i2,j2,2)+X(i2-1,j2,2))/2)/2)+1;
            c_1_4 = floor((X(i2-1,j2,3)+(X(i2,j2,3)+X(i2-1,j2,3))/2)/2)+1;
                
            if a2==237
                a2=236;
            end
            if b2==237
                b2=236;
            end
            if c2==237
                c2=236;
            end

            if a_medio==237
                a_medio=236;
            end
            if b_medio==237
                b_medio=236;
            end
            if c_medio==237
                c_medio=236;
            end
            
            if a_1_4==237
                a_1_4=236;
            end
            if b_1_4==237
                b_1_4=236;
            end
            if c_1_4==237
                c_1_4=236;
            end
            
            if a_3_4==237
                a_3_4=236;
            end
            if b_3_4==237
                b_3_4=236;
            end
            if c_3_4==237
                c_3_4=236;
            end
                
            if ((connected_matrix(a2,b2,c2)==0)&&(connected_matrix(a_medio,b_medio,c_medio)==0)&&(connected_matrix(a_1_4,b_1_4,c_1_4)==0)&&(connected_matrix(a_3_4,b_3_4,c_3_4)==0)&&(crossover==0))||((connected_matrix(a2,b2,c2)==0)&&(crossover==1))
                continue2 = 0;
            elseif i3==10
                X(i2,j2,1) = X(i2-1,j2,1);
                X(i2,j2,2) = X(i2-1,j2,2);
                X(i2,j2,3) = X(i2-1,j2,3);
                continue2 = 0;
            else
                i3 = i3+1;
            end
        end
    end
end
    
% save('C:\Users\Cliente\Documents\tesi\trajectories\random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE','-v7.3')
% save('C:\Users\Cliente\Documents\tesi\trajectories\random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE')
save('/Users/tiagocardoso/Documents/FISICA/Files_Matlab/tesi/trajectories/random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE','-v7.3')
% save('/Users/tiagocardoso/Documents/FISICA/Files_Matlab/tesi/trajectories/random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE')

status='random walk done'
    
end