function random_walk_function_little_spheres(TE)
% same as Random_walk but it's a function, unique file and for the little
% spheres

filename = 'spheres_500_54_256_I.mat';

num_spins = 10000;
% num_spins = 1000;
% num_spins = 3200;

% step_meter = (10^(-6))*18/50;
% step_meter = (10^(-6))*(3/50);
step_meter = (10^(-6))*(3/5);
% step_meter = (10^(-6))*(9/5);

Diam_meter = (30)*(10^(-6));

D = (2.3)*(10)^(-9);

% TE = (84.4)*(10)^(-3);
% TE = (74.4)*(10)^(-3);
% TE = (44.4)*(10)^(-3);
% TE = (3304.4)*(10)^(-3);
% TE = (0.002)*(10)^(-3);

load(filename, 'connected_matrix', 'Diam')

RES = length(connected_matrix);
dt = ((step_meter)^2)/(6*D);
numsteps = floor(TE/dt);

convert_coord_meter = Diam/Diam_meter;
% convert_coord_meter = 4040000;
step = step_meter*convert_coord_meter;

X = zeros(numsteps + 1,num_spins,3);
% X0 = 216*rand(num_spins,3)+20;
% X0 = 246*rand(num_spins,3)+5;
X0 = 56*rand(num_spins,3)+100;


for i1 = 1:num_spins
    continue1 = 1;
    while continue1
        a1 = floor(X0(i1,1))+1;
        b1 = floor(X0(i1,2))+1;
        c1 = floor(X0(i1,3))+1;
        
%         if a1==237
%             a1=236;
%         end
        
%         if a1==252
%             a1=251;
%         end
%         if b1==252
%             b1=251;
%         end
%         if c1==252
%             c1=251;
%         end
        
        if connected_matrix(a1,b1,c1) > 1000
            continue1 = 0;
        else
            X0(i1,1) = 56*rand+100;
            X0(i1,2) = 56*rand+100;
            X0(i1,3) = 56*rand+100;
%             X0(i1,1) = 246*rand+5;
%             X0(i1,1) = 56*rand+100;
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
            
%             if X(i2,j2,1) < 5
%                 crossover = 1;
%                 X(i2,j2,1) = X(i2,j2,1) + (fix((5-X(i2,j2,1))/246)+1)*246;
%             elseif X(i2,j2,1) >= 251
%                 crossover = 1;
%                 X(i2,j2,1) = X(i2,j2,1) - fix((X(i2,j2,1)-5)/246)*246;
%             end
            
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
            
            a1 = floor(X(i2-1,j2,1))+1;
            b1 = floor(X(i2-1,j2,2))+1;
            c1 = floor(X(i2-1,j2,3))+1;
            
            
%             if a2==252
%                 a2=251;
%             end
            if a2==237
                a2=236;
            end
            if b2==237
                b2=236;
            end
            if c2==237
                c2=236;
            end

            if a1==237
                a1=236;
            end
            if b1==237
                b1=236;
            end
            if c1==237
                c1=236;
            end
            
             
            if ((connected_matrix(a2,b2,c2)==connected_matrix(a1,b1,c1))&&(crossover==0))||((connected_matrix(a2,b2,c2)>1000)&&(crossover==1))
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
% save('/Users/tiagocardoso/Documents/FISICA/Files_Matlab/tesi/trajectories/random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE','-v7.3')
save('/Users/tiagocardoso/Documents/FISICA/Files_Matlab/tesi/trajectories/random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE')

status='random walk done'
    
end