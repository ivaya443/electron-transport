% close all
clear

counter = 1;
h_bar = double(1.054e-34);
m_e = double(0.067*9.1e-31);
q = double(1.602e-19);

a = double(5e-9);
b = double(h_bar^2/(2*m_e*a^2*q));
d = 0;

Nx=11-2; % This is the # of internal sites of the system. Exclude left and righ slice!!
Ny=19;
extra_sites = 400;

res_length = 300;

l_left = (Ny+1)*a;
l_right = (Ny+extra_sites+1)*a;

nnn = 2;
m = [1:Ny];
m_prim = [1:(Ny+extra_sites)];

H = zeros(Nx*Ny,Nx*Ny);
H_sub = zeros(Ny,Ny);
B = zeros(Ny,Ny);
B_left = zeros(Nx*Ny-Ny,Ny);
B_right = zeros(Nx*Ny-Ny,Ny);
one_min = -1*eye(Ny,Ny);
one_diag = eye(Ny,Ny);
cor_left = zeros(2*Ny,Ny);
cor_right = zeros(2*Ny+extra_sites,Ny+extra_sites);
Q = zeros(Nx*Ny+4*Ny+extra_sites,1);

E = double(4e-3);
k = double(sqrt(2*m_e*E*q/h_bar^2-m.^2*pi^2/(l_left)^2));
k_prim = double(sqrt(2*m_e*E*q/h_bar^2-m_prim.^2*pi^2/(l_right)^2));

for i = 1:Ny
    
    H_sub(i,i) = E - 4*b;
    B(i,i) = b;
    
    if i < Ny
        
        H_sub(i,i+1) = b;
        H_sub(i+1,i) = b;
        
    end
    
end

B_left = [B; B_left];
B_right = [B_right; B];

for i = 0:(Nx-1)
    
    x_range = (i*Ny+1):(i+1)*Ny;
    
    H(x_range,x_range) = H_sub;
    
    if i < (Nx-1)
        
        H(x_range,x_range+Ny) = B;
        H(x_range+Ny,x_range) = B;
        
    end
    
end

H = [B_left,H,B_right];

% aaa = [one_min,zeros(Ny,Nx*Ny+Ny)];
% bbb = [one_min,one_diag,zeros(Ny,(size(H,2))/2)];

one_top = [one_min,zeros(Ny,Nx*Ny+Ny);
    one_min,one_diag,zeros(Ny,Nx*Ny)];
one_bot = [zeros(Ny,Nx*Ny+Ny),one_min;
    zeros(Ny,Nx*Ny),one_diag,one_min];

H = [one_top; H; one_bot; zeros(extra_sites,Nx*Ny+2*Ny)];

for y = 1:2*Ny
    
    for m = 1:Ny
        
        if y <= Ny
            
            cor_left(y,m) = sin((pi*m*y*a)/l_left);
            
        end
        
        if y > Ny
            
            cor_left(y,m) = 1i*k(m)*a*sin((pi*m*(y-Ny)*a)/l_left);
            
        end
        
    end
    
end

for y = 1:2*Ny
    
    for m = 1:(Ny+extra_sites)
        
        if y <= Ny
            
            cor_right(y,m) = sin((pi*m*(y+extra_sites/2)*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
        if y > Ny
            
            cor_right(y,m) = 1i*k_prim(m)*a*sin((pi*m*(y-Ny+extra_sites/2)*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
    end
    
end

% Below are the extra rows for the extra sites!!

for y = 1:extra_sites
    
    for m = 1:(Ny+extra_sites)
        
        if y <= extra_sites/2
            
            cor_right(y+2*Ny,m) = sin((pi*m*y*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
        if y > extra_sites/2
            
            cor_right(y+2*Ny,m) = sin((pi*m*(y+Ny)*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
    end
    
end

cor_left = [cor_left; zeros(size(H,1)-2*Ny,Ny)];
cor_right = [zeros(size(H,1)-size(cor_right,1),Ny+extra_sites); cor_right];

H = [cor_left, H, cor_right];

for iii = 1:2*Ny
    
    if iii <= Ny
        
        Q(iii,1) = -sin(nnn*pi*iii*a/l_left);
        
    end
    
    if iii > Ny
        
        Q(iii,1) = 1i*k(nnn)*a*sin(nnn*pi*(iii-Ny)*a/l_left);
        
    end

end
 
psi_vec = linsolve(H,Q);

psi = reshape(psi_vec((Ny+1):(Nx*Ny+4*Ny+extra_sites-Ny-extra_sites)),[Ny,Nx+2]);

C_coeff = psi_vec(Nx*Ny+4*Ny+extra_sites-Ny-extra_sites+1:end,1);

psi_right = zeros(Ny+extra_sites,res_length);

for y = 1:(Ny+extra_sites)
    
    for x = (Nx+1):(Nx+res_length+1)
        
        sum = 0;
        
        for m = 1:(Ny+extra_sites)
            
            sum = sum + C_coeff(m)*sin(pi*m*y*a/l_right)*exp(1i*k_prim(m)*x*a);
            
        end
        
        psi_right(y,(x-Nx-0)) = sum;
        
    end
    
end

psi = [zeros(extra_sites/2 + 1,Nx+2); psi; zeros(extra_sites/2 + 1,Nx+2)];
psi_right = [zeros(1,res_length+1); psi_right; zeros(1,res_length+1)];
psi_full = [psi(:,1:(Nx + 1)),psi_right];

[d_x,d_y]=gradient(psi_full,a);

current_x = (h_bar/(2*m_e*1i))*(conj(psi_full).*d_x - psi_full.*conj(d_x));
current_y = (h_bar/(2*m_e*1i))*(conj(psi_full).*d_y - psi_full.*conj(d_y));

f111 = figure;
f222 = figure;

x = [0*1e9:a/1e-9:(Nx+1+res_length)*a/1e-9];
y = [-(Ny+extra_sites+1)*(5/2)*a/5e-9:a/1e-9:(Ny+extra_sites+1)*(5/2)*a/5e-9];

figure(f111)

[X,Y] = meshgrid(x,y);
Z = (psi_full.*conj(psi_full));
% Z = (psi_full+conj(psi_full))./2;
% Z = (psi_full-conj(psi_full))./2i;
% subplot(211)

surf_plot_1 = surf(X,Y,Z);

view(2)

hold on

% z_max = max(max(get(surf_plot_1,'Zdata')))
z_max = .1;
x_lgth = a*1e9*(Nx+1);
y_lgth = a*1e9*(Ny+1)/2;
line(0:x_lgth,y_lgth*ones(1,y_lgth+1),z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)
line(x_lgth*ones(1,x_lgth+1),y_lgth:2*y_lgth,z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)

line(0:x_lgth,-y_lgth*ones(1,y_lgth+1),z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)
line(x_lgth*ones(1,x_lgth+1),-y_lgth:-1:-2*y_lgth,z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)

hold off

ax = gca;
ax.SortMethod;
ax.SortMethod = 'childorder';
% ax.SortMethod = 'depth';

xlabel('x (nm)')
ylabel('y (nm)')

zlabel('|\Psi|^2')
title(['Hybrid |\Psi|^2 for E = ' num2str(E*1e3) ' meV'])

xlim([0 200])
ylim([-100 100])

% % plot(y,real(U_R_pos(:,1)))

figure(f222);

[X,Y] = meshgrid([0*1e9:a/1e-9:(Nx+1+res_length)*a/1e-9],[-(Ny+extra_sites+1)*(5/2)*a/5e-9:a/1e-9:(Ny+extra_sites+1)*(5/2)*a/5e-9]);

surf_plot_2 = surf(X,Y,sqrt(current_x.^2 + current_y.^2));

view(2)

hold on

% z_max = max(max(get(surf_plot_1,'Zdata')))
z_max = .1;
x_lgth = a*1e9*(Nx+1);
y_lgth = a*1e9*(Ny+1)/2;
line(0:x_lgth,y_lgth*ones(1,y_lgth+1),z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)
line(x_lgth*ones(1,x_lgth+1),y_lgth:2*y_lgth,z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)
line(0:x_lgth,-y_lgth*ones(1,y_lgth+1),z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)
line(x_lgth*ones(1,x_lgth+1),-y_lgth:-1:-2*y_lgth,z_max*ones(1,y_lgth+1),'Color','red','LineWidth',2)

hold off

ax = gca;
ax.SortMethod;
ax.SortMethod = 'childorder';
% ax.SortMethod = 'depth';
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('|J|')
title(['Hybrid |J| for E = ' num2str(E*1e3) ' meV'])


% hold on
% 
% x = [(Nx+1)*a/1e-9:a/1e-9:(Nx+1+res_length)*a/1e-9];
% y = [-(Ny+extra_sites+1)*(5/2)*a/5e-9:a/1e-9:(Ny+extra_sites+1)*(5/2)*a/5e-9];
% [X,Y] = meshgrid(x,y);
% Z = (psi_right.*conj(psi_right));
% % Z = (psi+conj(psi))./2;
% % Z = (psi-conj(psi))./2i;
% 
% surf(X,Y,Z)
% 
% hold off

xlim([0 200])
ylim([-100 100])