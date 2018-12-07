function bin_matrix = axons_256( N,PHI,RES )

% Fibers' Radii in coordinates (Splenium: 256 fibers) 8.8 for 80!!! /  8.1
% for 70%!!!

% r1 = 7.8*[3.51,3.51,2.70,2.70,2.43,2.43,2.03,2.03,2.03,1.76,1.76,1.76,1.76,1.62,1.62,1.62,1.62,1.62]';
% r2 = 7.8*[3.51,3.51,2.70,2.70,2.43,2.43,2.03,2.03,2.03,1.76,1.76,1.76,1.76,1.62,1.62,1.62,1.62,1.62,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22]';
% r3 = 7.8*[3.51,3.51,2.70,2.70,2.43,2.43,2.03,2.03,2.03,1.76,1.76,1.76,1.76,1.62,1.62,1.62,1.62,1.62,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,...
%     1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,...
%     0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,...
%     0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,...
%     0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68]';
% r4 = 7.8*[3.51,3.51,2.70,2.70,2.43,2.43,2.03,2.03,2.03,1.76,1.76,1.76,1.76,1.62,1.62,1.62,1.62,1.62,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,1.49,...
%     1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,...
%     0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,...
%     0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,...
%     0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.54,0.54,0.54,0.54,0.54,0.54,...
%     0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,...
%     0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,...
%     0.41,0.41,0.41,0.41,0.41,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.14,0.14,0.14,0.14,0.14]';


% Fibers' Radii in coordinates (Posterior Body: 256 fibers)

r1 = 6.7*[4.05,4.05,3.38,3.38,3.24,3.24,3.24,2.97,2.97,2.84,2.84,2.84,2.84]';
r2 = 6.7*[4.05,4.05,3.38,3.38,3.24,3.24,3.24,2.97,2.97,2.84,2.84,2.84,2.84,2.70,2.70,2.57,2.57,2.57,2.30,2.30,2.30,2.16,2.16,2.03,2.03,1.89,1.89,1.89,1.89,1.76,1.76,1.76,1.76,1.76,1.76]';
r3 = 6.7*[4.05,4.05,3.38,3.38,3.24,3.24,3.24,2.97,2.97,2.84,2.84,2.84,2.84,2.70,2.70,2.57,2.57,2.57,2.30,2.30,2.30,2.16,2.16,2.03,2.03,1.89,1.89,1.89,1.89,1.76,1.76,1.76,1.76,1.76,1.76,...
    1.62,1.62,1.62,1.62,1.62,1.49,1.49,1.49,1.49,1.49,1.49,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,...
    1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95]';
r4 = 6.7*[4.05,4.05,3.38,3.38,3.24,3.24,3.24,2.97,2.97,2.84,2.84,2.84,2.84,2.70,2.70,2.57,2.57,2.57,2.30,2.30,2.30,2.16,2.16,2.03,2.03,1.89,1.89,1.89,1.89,1.76,1.76,1.76,1.76,1.76,1.76,...
    1.62,1.62,1.62,1.62,1.62,1.49,1.49,1.49,1.49,1.49,1.49,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,1.22,...
    1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,1.08,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,...
    0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.81,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,...
    0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.68,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,...
    0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.54,0.41,0.41,0.41,...
    0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,0.27,...
    0.27,0.27,0.27,0.14,0.14,0.14]';
    


function f = mono_sphere_potential( C )
% calculates the potential in the case of non-penetrating polycircles

S = size(C);
num_circles = S(1);
f = 0;

for i=1:num_circles
    for j=1:num_circles
        if i < j
% The best up to now? No
            f = f + (10000)*((1 - ((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)/((r4(i)+r4(j))^2))^4)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r4(i)+r4(j)));
            
% ?
%             f = f + (10000)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)))^4)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)));
% ?
%             f = f + (10000)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)))^2)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)));
% Original
%             f = f + (1/2)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)))^2)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)));
% The New Best:
%             f = f + (10000)*((1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)))^4)*heaviside(1 - (((C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2)^(1/2))/(r(i)+r(j)))^4;
        end
    end
end

end

% Initialization of the centers' coordinates
x0 = (RES-2*r1).*rand(length(r1),1)+r1;
y0 = (RES-2*r1).*rand(length(r1),1)+r1;
X0 = x0;
X0(:,2) = y0;

% Bounds
lbx = r1;
lby = r1;
lb = lbx;
lb(:,2) = lby;

ubx = RES*ones(length(r1),1)-r1;
uby = RES*ones(length(r1),1)-r1;
ub = ubx;
ub(:,2) = uby;

% Unused constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Minimization Procedure (= final coordinates of the spheres' centers)
X0 = fmincon(@mono_sphere_potential,X0,A,b,Aeq,beq,lb,ub);

opt_run='run 1 complete'

% 2nd run of optimization algorithm

base = length(r1);
dimen = length(r2)-length(r1);
x0 = zeros(dimen,1);
y0 = zeros(dimen,1);

for p1=1:dimen
    continue21=1;
    while continue21
        x0(p1) = (RES-2*r4(base+p1))*rand + r4(base+p1);
        y0(p1) = (RES-2*r4(base+p1))*rand + r4(base+p1);
        
        for q1=1:base
            continue21=0;
            if (X0(q1,1)-x0(p1))^2 + (X0(q1,2)-y0(p1))^2 < (r4(q1)+r4(base+p1))^2
                continue21=1;
                break
            end
        end
    end
end

x0(:,2)=y0;
X0 = vertcat(X0,x0);

lbx = r2;
lby = r2;
lb = lbx;
lb(:,2) = lby;

ubx = RES*ones(length(r2),1)-r2;
uby = RES*ones(length(r2),1)-r2;
ub = ubx;
ub(:,2) = uby;

X0 = fmincon(@mono_sphere_potential,X0,A,b,Aeq,beq,lb,ub);

opt_run='run 2 complete'

% 3rd run of optimization algorithm

base = length(r2);
dimen = length(r3)-length(r2);
x0 = zeros(dimen,1);
y0 = zeros(dimen,1);

for p1=1:dimen
    continue21=1;
    while continue21
        x0(p1) = (RES-2*r4(base+p1))*rand + r4(base+p1);
        y0(p1) = (RES-2*r4(base+p1))*rand + r4(base+p1);
        
        for q1=1:base
            continue21=0;
            if (X0(q1,1)-x0(p1))^2 + (X0(q1,2)-y0(p1))^2 < (r4(q1)+r4(base+p1))^2
                continue21=1;
                break
            end
        end
    end
end

x0(:,2)=y0;
X0 = vertcat(X0,x0);

lbx = r3;
lby = r3;
lb = lbx;
lb(:,2) = lby;

ubx = RES*ones(length(r3),1)-r3;
uby = RES*ones(length(r3),1)-r3;
ub = ubx;
ub(:,2) = uby;

X0 = fmincon(@mono_sphere_potential,X0,A,b,Aeq,beq,lb,ub);

opt_run='run 3 complete'

% 4th run of optimization algorithm

base = length(r3);
dimen = length(r4)-length(r3);
x0 = zeros(dimen,1);
y0 = zeros(dimen,1);

for p1=1:dimen
    continue21=1;
    while continue21
        x0(p1) = (RES-2*r4(base+p1))*rand + r4(base+p1);
        y0(p1) = (RES-2*r4(base+p1))*rand + r4(base+p1);
        
        for q1=1:base
            continue21=0;
            if (X0(q1,1)-x0(p1))^2 + (X0(q1,2)-y0(p1))^2 <= (r4(q1))^2
                continue21=1;
                break
            end
        end
    end
end

x0(:,2)=y0;
X0 = vertcat(X0,x0);

lbx = r4;
lby = r4;
lb = lbx;
lb(:,2) = lby;

ubx = RES*ones(length(r4),1)-r4;
uby = RES*ones(length(r4),1)-r4;
ub = ubx;
ub(:,2) = uby;

X = fmincon(@mono_sphere_potential,X0,A,b,Aeq,beq,lb,ub);

opt_run='run 4 complete'

% Binarization of the volumes
x=0.5:1:(RES-0.5);
y=x;
bin_matrix=zeros(RES,RES,RES);

for n=1:N
    for k=1:RES
        for l=1:RES
            if ((x(k)-X(n,1))^2+(y(l)-X(n,2))^2 > (0.74*(r4(n)))^2)&&((x(k)-X(n,1))^2+(y(l)-X(n,2))^2 <= ((r4(n)))^2)
                bin_matrix(k,l,:)=1;
            end
            if (x(k)-X(n,1))^2+(y(l)-X(n,2))^2 <= (0.74*(r4(n)))^2
                bin_matrix(k,l,:)=2;
            end
        end
    end
end

function outvec=demyelin(invec)
    
    outvec = invec;
    [numrow, numcol, numfloor] = size(invec);
    
    for l1=1:numrow
        for l2=1:numcol
            
            a1=invec(l1,l2,1);
            b1=invec(l1,l2,2);
            c1=invec(l1,l2,3);
    
            if a1==1
                nr=[a1,a1+1];
            elseif a1==RES
                nr=[a1-1,a1];
            else
                nr=[a1-1,a1,a1+1];
            end
            
            if b1==1
                nc=[b1,b1+1];
            elseif b1==RES
                nc=[b1-1,b1];
            else
                nc=[b1-1,b1,b1+1];
            end
            
%             if c1==1
%                 nh=[c1,c1+1];
%             elseif c1==RES
%                 nh=[c1-1,c1];
%             else
%                 nh=[c1-1,c1,c1+1];
%             end
            
%             acrescimo para restringir a inflamacao no sentido
%             longitudinal
            
%             if c1<(RES/2)
%                 nh=[c1,c1+1];
%             else
%                 nh=[c1-1,c1];
%             end
            
            if (a1<85)&(b1>85)&(b1<170)
                if c1<176
                    nh=[c1,c1+1];
                else
                    nh=[c1-1,c1];
                end
            elseif (a1>85)&(a1<170)&(b1<85)
                if c1<80
                    nh=[c1,c1+1];
                else
                    nh=[c1-1,c1];
                end
            elseif (a1>85)&(a1<170)&(b1>170)
                if c1<100
                    nh=[c1,c1+1];
                else
                    nh=[c1-1,c1];
                end
            elseif (a1>170)&(b1>85)&(b1<170)
                if c1<156
                    nh=[c1,c1+1];
                else
                    nh=[c1-1,c1];
                end
            elseif (a1>85)&(a1<170)&(b1>85)&(b1<170)
                if c1<128
                    nh=[c1,c1+1];
                else
                    nh=[c1-1,c1];
                end
            else
                if c1==1
                    nh=[c1,c1+1];
                elseif c1==RES
                    nh=[c1-1,c1];
                else
                    nh=[c1-1,c1,c1+1];
                end
                
            end    
            
            continue1 = 1;
            count = 0;
            
            while continue1 && (count<1)
                
                a_vic = datasample(nr,1);
                b_vic = datasample(nc,1);
                c_vic = datasample(nh,1);
                
                if (bin_matrix(a_vic,b_vic,c_vic)==1)&&((a_vic==a1)+(b_vic==b1)+(c_vic==c1)~=3)
                    outvec(l1,l2,:) = [a_vic,b_vic,c_vic];
                    continue1 = 0; 
                else
                    count = count+1;
                end
            end
        end
           
    end
        
end

function labelling(ctl)
    
    [numrows, numcols, numfloors] = size(ctl);
    
    for k1=1:numrows
        for k2=1:numcols
            bin_matrix(ctl(k1,k2,1),ctl(k1,k2,2),ctl(k1,k2,3))=0;
        end
    end
    
end


if PHI ~= 1
    
%     Nr of focal pts of demyelination per axon (300 for 30%, 600 for 60%)
%     num_focal = 2400;
    num_focal = 300;
    
    chosen = zeros(num_focal,N,3);
    
    for ordo=1:num_focal
        
        cosines = 2*rand(N,1)-1;
        
        pos_sign = randi(2,N,1);
        the_sign = [-1,1]';
        sines = ((1-(cosines).^2).^(1/2)).*the_sign(pos_sign);
    
        height = randi(256,N,1);
        
        POSx = X(:,1)+r4.*cosines;
        POSy = X(:,2)+r4.*sines;
        
        CELLx = floor(POSx)+1;
        CELLy = floor(POSy)+1;
        
        if CELLx > RES
            CELLx = RES;
        end
        if CELLy > RES
            CELLy = RES;
        end
    
        for n1=1:N
            
            chosen(ordo,n1,:) = [CELLx(n1),CELLy(n1),height(n1)];
        
            if (cosines(n1)>=0)&&(sines(n1)>=0)&&(bin_matrix(CELLx(n1),CELLy(n1),height(n1))~=1)
                chosen(ordo,n1,:) = [CELLx(n1)-1,CELLy(n1)-1,height(n1)];
            end
            if (cosines(n1)>=0)&&(sines(n1)<0)&&(bin_matrix(CELLx(n1),CELLy(n1),height(n1))~=1)
                chosen(ordo,n1,:) = [CELLx(n1)-1,CELLy(n1)+1,height(n1)];
            end
            if (cosines(n1)<0)&&(sines(n1)>=0)&&(bin_matrix(CELLx(n1),CELLy(n1),height(n1))~=1)
                chosen(ordo,n1,:) = [CELLx(n1)+1,CELLy(n1)-1,height(n1)];
            end
            if (cosines(n1)<0)&&(sines(n1)<0)&&(bin_matrix(CELLx(n1),CELLy(n1),height(n1))~=1)
                chosen(ordo,n1,:) = [CELLx(n1)+1,CELLy(n1)+1,height(n1)];
            end
        end
    end
    
    fila = chosen;
    M = sum(bin_matrix(:)==1);
    m = M;
    proportion =m/M;
    
    while proportion > PHI
        
        labelling(fila);
        proportion = (sum(bin_matrix(:)==1))/M
        fila = demyelin(fila);
        
    end
    
    
    
    
end



end