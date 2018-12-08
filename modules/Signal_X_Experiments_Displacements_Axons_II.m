%{ 

Main Program for the Signal Calculation Block (axons) including: 

(1) multiple experiments and (2) signal calculation based on the displacements from initial position. 
Made for the Axons simulations. Please note that experimental time ("delta") is fixed. Consequently, the values of variable 
bfield depend exclusively on the gradient values employed.

(C) 2018 Tiago Monteiro Cardoso
Last update: 2018-12-08

%}


clear

% parpool

num_experiments = 15;
% num_experiments = 1;

file_geometry = 'interaxon_post256_PHI_80_DEM_60_XIII.mat';

% direction of application of gradient
grad_direction = [0,1,0];
deltachi = -(7.8)*(10^(-8));
B0 = 9.4;
% B0 = 1;
gamma = 267500000;
pulse = (4.4)*(10)^(-3);
delta = (80)*(10)^(-3);
% delta = (40)*(10)^(-3);

data = zeros(num_experiments,3);

% bfield = (10^6)*[100,500,1000,1500,2000,2500,3000,3500,4000,5000,6000,7000,8000,9000,10000];
bfield = (10^6)*[100,500,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000];


for i4=1:num_experiments
    
    random_walk_function_axons_II
    
%     load('C:\Users\Cliente\Documents\tesi\trajectories\random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE')
    load('/Users/tiagocardoso/Documents/FISICA/Files_Matlab/tesi/trajectories/random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE')
        
    [numsteps,N,coordinates] = size(X);
        
    Xsel = X(:,:,1);
    Ysel = X(:,:,2);
    Zsel = X(:,:,3);
  
    clear X

%     LENGTH PDF ROUTINE
%     sizeZ = size(Zsel);
%     
%     if i4==1
%         Ztrajectories = zeros(numsteps,max_spins*num_experiments);
%         trajcount=1;
%     end
%     
%     Ztrajectories(:,trajcount:trajcount+sizeZ(2)-1) = Zsel;
%     trajcount = trajcount+sizeZ(2);
%     LENGTH PDF ROUTINE
    
    steps_pulse = floor(pulse/dt);
    size_time = size(Xsel);
    timeline = zeros(size_time(1)-1,size_time(2));
    
    for timeindex=1:size_time(1)-1
        
        if timeindex <= steps_pulse
            timeline(timeindex,:) = 1;
        end
        if timeindex >= size_time(1)-steps_pulse
            timeline(timeindex,:) = -1;
        end
    end
    
    gradline = ones(size(Xsel));
    hec = floor(numsteps/2)+1;
    gradline(hec:numsteps,:)=-1;
    
    ghat = grad_direction/norm(grad_direction);
    k = (bfield(i4)/((gamma^2)*(pulse^2)*(delta-pulse/3)))^(1/2);
    gx = k*ghat(1);
    gy = k*ghat(2);
    gz = k*ghat(3);
    
    Bbase = B0;
    
    deltaphiwithout = sum(gamma*dt*((gx/convert_coord_meter)*cumsum(diff(Xsel)).*timeline+(gy/convert_coord_meter)*cumsum(diff(Ysel)).*timeline+(gz/convert_coord_meter)*cumsum(diff(Zsel)).*timeline));
%     deltaphiwith = deltaphiwithout + sum(gamma*dt*internal_gradient(Xsel,Ysel,Zsel,file_geometry,deltachi,Bbase).*gradline);
    
%     E1 = abs(sum(exp(1i*deltaphiwithout))/i2);
%     E1 = sum(real(exp(1i*deltaphiwithout)))/i2;
    E1 = sum(cos(deltaphiwithout))/N;

%     signal intensity:
    E2 = abs(sum(exp(1i*deltaphiwithout))/N);

%     E2 = sum(cos(deltaphiwith))/N;
    
    data(i4,1) = bfield(i4);
    data(i4,2) = E1;
    data(i4,3) = E2;
    
    Experiment_Done = i4
    
end

% Ztrajectories(:,all(Ztrajectories==0))=[];
% (LENGTH PDF ROUTINE^)

x1=data(:,1)/1000000000;
y1=data(:,2);
y2=data(:,3);

% fun = @(x,x1)mlf(x(1)/1000,1,-(x(2).*x1).^x(3),11);
fun = @(x,x1)mlf(x(1),1,-(x(2).*x1).^x(3),11);

x0=[0.5 0.1 0.5];
% x0=[0.5 0.01 0.5];

% x0=[45 0.01 0.5];

lb=[];
ub=[];
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

x = lsqcurvefit(fun,x0,x1,y1,lb,ub,options)

h(1) = figure;
plot(x1,y1,'ko')
hold on
t=linspace(x1(1),x1(end),100);
plot(t,fun(x,t),'k-')
xlabel('b-value (1000 s/mm^{2})')
ylabel ('S/S0')
legend ('Experimental points' , 'Fit' )
% set(gca, 'YScale', 'log')
title(['Signal (without int. gradient) - \alpha = ' num2str(x(1),'%.2f') ' / D = ' num2str(x(2),'%.2f') ' x 10^{-9} m^{2}.s^{-1} / \gamma = ' num2str(x(3),'%.2f')])

x = lsqcurvefit(fun,x0,x1,y2,lb,ub,options)

h(2) = figure;
plot(x1,y2,'ko')
hold on
t=linspace(x1(1),x1(end),100);
plot(t,fun(x,t),'k-')
xlabel('b-value (1000 s/mm^{2})')
ylabel ('S/S0')
legend ('Experimental points' , 'Fit' )
% set(gca, 'YScale', 'log')
title(['Signal (absolute value) - \alpha = ' num2str(x(1),'%.2f') ' / D = ' num2str(x(2),'%.2f') ' x 10^{-9} m^{2}.s^{-1} / \gamma = ' num2str(x(3),'%.2f')])



% x1=data(:,1)/1000000;
% y1=log(data(:,2));
% y2=log(data(:,3));

% h(1) = figure;
% plot(x1,y1,'ko')
% hold on
% xlabel('b-value (s/mm^{2})')
% ylabel ('Ln(S/S0)')
% fun = @(x,x1)-(x(1).*x1).^x(2);
% x0=[1 1];
% % fun = @(x,x1)exp(-(x(1).*x1).^x(2));
% % x0=[0.5 0.5];
% lb=[];
% ub=[];
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
% x = lsqcurvefit(fun,x0,x1,y1,lb,ub,options);
% % Diff_without = x(1)*1000;
% Diff_without = x(1)*10000;
% t=linspace(x1(1),x1(end),100);
% plot(t,fun(x,t),'k-')
% legend ('Experimental points' , 'Fit' )
% title(['Signal (without int. gradient) - D = ' num2str(Diff_without,'%.2f') ' x 10^{-10} / \gamma = ' num2str(x(2),'%.2f')])


% h(2) = figure;
% plot(x1,y2,'ro')
% hold on
% xlabel ( 'b-value (s/mm^{2})' )
% ylabel ('Ln(S/S0)')
% fun = @(x,x1)-(x(1).*x1).^x(2);
% x0=[1 1];
% % fun = @(x,x1)x(1)-(x(2).*x1).^x(3);
% % x0=[1 1 1];
% % fun = @(x,x1)exp(-(x(1).*x1).^x(2));
% % x0=[0.5 0.5];
% lb=[];
% ub=[];
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
% x = lsqcurvefit(fun,x0,x1,y2,lb,ub,options);
% Diff_with = x(1)*10000;
% t=linspace(x1(1),x1(end),100);
% plot(t,fun(x,t),'r-')
% legend ('Experimental points' , 'Fit' )
% title(['Signal (with int. gradient) - D = ' num2str(Diff_with,'%.2f') ' x 10^{-10} / \gamma = ' num2str(x(2),'%.2f')])
% 
% 
savefig(h,'Signal two figures.fig')
% close(h)
% 
% % q=sqrt(bfield'/(delta-pulse/3))/2*pi

save('spheres_results.mat','data','B0','grad_direction','pulse','deltachi','N')
