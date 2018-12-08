%{ 

Main Program for the Signal Calculation Block (spheres) including: 

(1) multiple experiments and (2) signal calculation based on the displacements from initial position. 
Made for the spheres simulations. Please note that experimental time (variable "delta") is NOT fixed. Consequently, 
the values of variable bfield depend on both experimental time AND the gradient values employed.

(C) 2018 Tiago Monteiro Cardoso
Last update: 2018-12-08

%}

clear

% parpool

% num_experiments = 24;
num_experiments = 48;

% num_experiments = 16;

file_geometry = 'spheres_500_70_256_I.mat';

% direction of application of gradient
grad_direction = [1,0,0];
deltachi = (1)*(10^(-7));
% deltachi = 0.0000003;
B0 = 9.4;

gamma = 267500000;
pulse = (4.4)*(10)^(-3);
% pulse = (1)*(10)^(-3);
% delta = (80)*(10)^(-3);
% delta = (3300)*(10)^(-3);
% delta = (40)*(10)^(-3);
% delta = (0.001)*(10)^(-3);

% data = zeros(num_experiments,3);
data = zeros(num_experiments,4);

% bfield = (10^6)*[100,500,1000,1500,2000,2500,3000,3500,4000,5000,6000,7000,8000,9000,10000];
% bfield = (10^6)*[0,100,200,400,600,800,1000,1250,1500,1873,2200,2500,3000,4000,5000,6500,7492,8500,10000,12000,14000,16000,16857,17800];
% bfield = (10^6)*[100,500,1000,1500,2000,2500,3000,3500,4000,5000,6000,8000,10000,15000,20000,25000];

bfield = (10^6)*[100,103,106,108,110,113,116,118,300,310,320,330,340,350,360,370,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,25000];

% delta = (10^(-3))*[5.6,30,40,60,70,80,100,120];
delta = (10^(-3))*[30,40,50,60,70,80,100,120];


for i4=1:num_experiments
    
    random_walk_function_little_spheres(pulse+delta(mod(i4-1,8)+1))
%     random_walk_function_little_spheres
    
    load('/Users/tiagocardoso/Documents/FISICA/Files_Matlab/tesi/trajectories/random_walk_unique.mat', 'X','convert_coord_meter','RES','dt','TE')
        
    [numsteps,N,coordinates] = size(X);
        
    Xsel = X(:,:,1);
    Ysel = X(:,:,2);
    Zsel = X(:,:,3);
    
%     Routine for the computation of mean squared displacement
    
    A = X(:,:,1);
    Differences = diff(A);
    Positive_crossover = sum(Differences(:,:)>100);
    Negative_crossover = sum(Differences(:,:)<-100);
    Total_crossover = Positive_crossover - Negative_crossover;
    squared_displacement = ((A(numsteps,:)-A(1,:)-216*Total_crossover)/convert_coord_meter).^2;
    MSD_of_experiment = mean(squared_displacement);
    
%     End of Routine
    
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
%     k = (bfield(i4)/((gamma^2)*(pulse^2)*(delta-pulse/3)))^(1/2);
    k = (bfield(i4)/((gamma^2)*(pulse^2)*(delta(mod(i4-1,8)+1)-pulse/3)))^(1/2);
    gx = k*ghat(1);
    gy = k*ghat(2);
    gz = k*ghat(3);
    
    Bbase = B0;
    
    deltaphiwithout = sum(gamma*dt*((gx/convert_coord_meter)*cumsum(diff(Xsel)).*timeline+(gy/convert_coord_meter)*cumsum(diff(Ysel)).*timeline+(gz/convert_coord_meter)*cumsum(diff(Zsel)).*timeline));
%     deltaphiwith = deltaphiwithout + sum(gamma*dt*internal_gradient(Xsel,Ysel,Zsel,file_geometry,deltachi,Bbase).*gradline);
    
%     E1 = abs(sum(cos(deltaphiwithout))/N);
%     E1 = abs(sum(exp(1i*deltaphiwithout))/N);
%     E1 = sum(real(exp(1i*deltaphiwithout)))/i2;
    E1 = sum(cos(deltaphiwithout))/N;
    

%     E2 = sum(cos(deltaphiwith))/N;
    
%     data(i4,1) = bfield(i4);
    data(i4,1) = bfield(i4)/(delta(mod(i4-1,8)+1)-pulse/3);
    
    data(i4,2) = E1;
%     data(i4,3) = E2;
    data(i4,3) = delta(mod(i4-1,8)+1)-pulse/3;
    data(i4,4) = MSD_of_experiment;
    
    Experiment_Done = i4
    
end

% Ztrajectories(:,all(Ztrajectories==0))=[];
% (LENGTH PDF ROUTINE^)

% x1=data(:,1)/1000000000;
% y1=data(:,2);
% y2=data(:,3);
% MSD=mean(data(:,4))

x1(:,1)=data(:,1)/1000000000;
y1=data(:,2);
x1(:,2)=data(:,3);
MSD=mean(data(:,4))

% q=(sqrt(data(:,1)/(delta-pulse/3))/(2*pi))*(5*10^(-5));
% delta_bar=(delta-pulse/3);

% fun = @(x,x1)mlf(x(1),1,-((x(2).*x1).^(2*x(3))).*delta_bar.^x(1),11);
% fun = @(x,x1)mlf(x(1),1,-(x(2).*x1),11);
% fun = @(x,x1)mlf(x(1),1,-(x(2).*x1).^x(3),11);
fun = @(x,x1)mlf(x(1),1,-x(2).*(x1(:,1).^x(3)).*(x1(:,2).^x(1)),11);

x0=[0.5 1.5 0.5];
% x0=[0.9 0.5 0.9];
% x0=[0.5 1.5];

lb=[];
ub=[];
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

x = lsqcurvefit(fun,x0,x1,y1,lb,ub,options)

figure
x = 0:20:800;
y = 0.02:0.005:0.125;
[X,Y] = meshgrid(x,y);
Z = mlf(0.86,1,-0.92.*(X.^0.93).*(Y.^0.86),11);
s=surf(X,Y,Z);
s.FaceAlpha=0.5;
s.FaceLighting='gouraud';
hold on
p=plot3(x1(:,1),x1(:,2),y1,'ko');
p.MarkerFaceColor='black';


figure
plot3(x1(:,1),x1(:,2),y1,'ko')
% hold on
% s=linspace(x1(1,1),x1(36,1),100);
% t=linspace(x1(1,2),x1(36,2),100);
% plot3(t,fun(x,t),'k-')
% 
% h(1) = figure;
% plot(x1,y1,'ko')
% hold on
% t=linspace(x1(1),x1(end),100);
% plot(t,fun(x,t),'k-')
% xlabel('b-value (1000 s/mm^{2})')
% ylabel ('S/S0')
% legend ('Experimental points' , 'Fit' )
% % set(gca, 'YScale', 'log')
% title(['Signal (without int. gradient) - \alpha = ' num2str(x(1),'%.2f') ' / D = ' num2str(x(2),'%.2f') ' x 10^{-9} m^{2}.s^{-1} / \gamma = ' num2str(x(3),'%.2f')])

% fun = @(x,x1)mlf(x(1),1,-(x(2).*x1).^x(3),11);
% x = lsqcurvefit(fun,x0,x1,y2,lb,ub,options)
% 
% h(2) = figure;
% plot(x1,y2,'ro')
% hold on
% plot(t,fun(x,t),'r-')
% xlabel('b-value (1000 s/mm^{2})')
% ylabel ('S/S0')
% legend ('Experimental points' , 'Fit' )
% % set(gca, 'YScale', 'log')
% title(['Signal (with int. gradient) - \alpha = ' num2str(x(1),'%.2f') ' / D = ' num2str(x(2),'%.2f') ' x 10^{-9} m^{2}.s^{-1} / \gamma = ' num2str(x(3),'%.2f')])
% 
% savefig(h,'Signal two figures ML.fig')
% close(h)
% 
save('spheres_results_ML.mat','data','B0','grad_direction','pulse','deltachi','N')
