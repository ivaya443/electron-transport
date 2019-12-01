% close all
clear

%% Constants:

counter = 1; % Old bit of code!
h_bar = double(1.054e-34);
m_e = double(0.067*9.1e-31);
q = double(1.602e-19);

%% Geometries:

a = double(5e-9);
b = double(h_bar^2/(2*m_e*a^2*q));
% d = 0;

Nx = 11-2; % This is # of internal sites in the system. Exclude left and right slice!
Ny = 19; % This is the left lead (smaller wire) width!!
extra_sites = 400; % Multiple of 2!!

res_length = 30; % Length of reservoir!!

l_left = (Ny + 1)*a; % Left lead!
l_right = (extra_sites + Ny+2*(Nx-1) + 1)*a; % Right lead; think how the triangular form is created!

%% Modes:

nnn = 2;
m = 1:Ny; % Mode # must match sites in y-direction!
m_prim = 1:(2*(Nx-1)+Ny+extra_sites); % Check carefully!

%% Energy & wave numbers:

E = double(4e-3);
k = double(sqrt(2*m_e*E*q/h_bar^2-m.^2*pi^2/(l_left)^2));
k_prim = double(sqrt(2*m_e*E*q/h_bar^2-m_prim.^2*pi^2/(l_right)^2));

%% Initialization:

% H = zeros(3*Ny+4*Nx+2*Nx+extra_sites,3*Ny+4*Nx+2*Nx+extra_sites); % Index
% wrong!!
H_sub = zeros(Ny,Ny); % Atom: smallest internal points slice!
B = zeros(Ny,Ny); % Transfer integral between slices!
B_left = zeros(Ny,Ny); % Transfer integral for left lead!
% B_right = zeros(2*Nx,2*Nx); % -||- for right lead!
one_min = -1*eye(Ny,Ny); % Minus identity!
one_diag = eye(Ny,Ny); % Identity!
cor_left = zeros(2*Ny,Ny); % Left lead!
cor_right = zeros(2*(Ny+2*(Nx-1))+extra_sites,2*(Nx-1)+Ny+extra_sites); % Right lead!
% Q = zeros(3*Ny+4*Nx+2*Nx+extra_sites,1); % Right-hand side of Schrödinger equation! Check index!

%% Left lead:

for y = 1:2*Ny
    
    for m = 1:Ny
        
        if y <= Ny % Upper half, i.e. wave function matching!
            
            cor_left(y,m) = sin((pi*m*y*a)/l_left);
            
        end
        
        if y > Ny % Bottom half, i.e. derivatives matching! Check sines!
            
            cor_left(y,m) = 1i*k(m)*a*sin((pi*m*(y-Ny)*a)/l_left);
            
        end
        
    end
    
end

for i = 1:Ny
    
    H_sub(i,i) = E - 4*b;
    B(i,i) = b;
    
    if i < Ny
        
        H_sub(i,i+1) = b;
        H_sub(i+1,i) = b;
        
    end
    
end

B_left = [B_left,B]; % Zeros and a B connection with left lead!

lead_l = vertcat([cor_left,[one_min, zeros(Ny,Ny);
                            one_min, one_diag]],[B_left, H_sub]);

% lead_l = horzcat([cor_left; zeros(Ny,Ny)],[one_min, zeros(Ny,Ny);
%                                            one_min, one_diag;
%                                            B_left, H_sub]);

%% Internal sites:

H = lead_l;

for i = 1:(Nx-1)
    
%     tmp = repmat({B},i,1); % Repeat B blocks
%     B_tmp = [zeros(2*i,1) , blkdiag(tmp{:}), zeros(2*i,1)];
    
    B_tmp = [zeros(Ny+2*(i-1),1) , b*eye(Ny+2*(i-1)), zeros(Ny+2*(i-1),1)]; % Build custom B blocks
    
%     tmp = repmat({H_sub},2*i,1); % Repeat H_sub blocks
%     H_tmp = blkdiag(tmp{:});

    H_tmp = full(gallery('tridiag',Ny+2*i,b,(E - 4*b),b)); % Build custom H_sub blocks
    
    zero_block = zeros(size(H,1) - (Ny + 2*(i-1)), size(H_tmp,1));
    
    H = [H, [zero_block; B_tmp]; % Should be possible to pre-allocate for speedup!
         zero_block', B_tmp', H_tmp];
    
end

%% Right lead:

opening = Ny+2*(Nx-1); % How large is the right lead (opening to reservoir)?!
last_slice = eye(opening); % Final slice block size!

for y = 1:2*opening
    
    for m = 1:(2*(Nx-1)+Ny+extra_sites)
        
        if y <= opening % Upper half, just the wave functions matching!
            
            cor_right(y,m) = sin((pi*m*(y+extra_sites/2)*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
        if y > opening % Bottom half (excluding extra points) the derivatives matching!
            
            cor_right(y,m) = 1i*k_prim(m)*a*sin((pi*m*(y-opening+extra_sites/2)*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
    end
    
end

for y = 1:extra_sites % Last rows in cor_right are for the extra points!
    
    for m = 1:(2*(Nx-1)+Ny+extra_sites)
        
        if y <= extra_sites/2 % Extra sites closer to origo!
            
            cor_right(y+2*opening,m) = sin((pi*m*y*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
        if y > extra_sites/2 % Extra sites farther away from origo!
            
            cor_right(y+2*opening,m) = sin((pi*m*(y+opening)*a)/l_right)*exp(1i*k_prim(m)*(Nx+1)*a);
            
        end
        
    end
    
end

%% Full matrix H & right-hand side vector Q:

B_right = b*last_slice;
one_min_right = -1*last_slice;
H_size = size(H,1);
B_r_size = size(B_right,1);
BBB = [zeros(B_r_size), one_min_right;
       last_slice, one_min_right;
       zeros(extra_sites,2*B_r_size)];

H = [H, [zeros(H_size-B_r_size,B_r_size); B_right], zeros(H_size,size(cor_right,2));
     zeros(size(cor_right,1),H_size-B_r_size), BBB, cor_right];
 
Q = zeros(size(H,1),1);

for iii = 1:2*Ny
    
    if iii <= Ny
        
        Q(iii,1) = -sin(nnn*pi*iii*a/l_left);
        
    end
    
    if iii > Ny
        
        Q(iii,1) = 1i*k(nnn)*a*sin(nnn*pi*(iii-Ny)*a/l_left);
        
    end
    
end

%% Solver:

% psi = reshape(psi_vec(Ny+1),[Ny,Nx+2]);

% p = [ 1 1 2 2 3 3 3 3 4 4 4 4];
% aq = [0, 0, 1, 1;
%       1, 1, 1, 1;
%       1, 1, 1, 1;
%       0, 0, 1, 1];
% bam = tril(ones(4)).'
% aq(aq > 0) = p;

% bq = [[zeros((opening-Ny),1);ones(Ny,1);zeros((opening-Ny),1)], [tri_up; tri_dn], ones(opening,1)];
% bq = [[zeros, tri_up, ];
%       ones(Ny,Nx+2)]];

% psi_vec = linsolve(H,Q)';
psi_vecc = linsolve(H,Q); % Important correction, have psi_vec standing up!!

tri_up = fliplr(tril(ones((Nx-1))));
tri_dn = flipud(tri_up);

system = [zeros(size(tri_up,1),2), tri_up, ones(size(tri_up,1),1);
          ones(Ny,Nx+2);
          zeros(size(tri_dn,1),2), tri_dn, ones(size(tri_dn,1),1)];

% % system(system==0) = -1;
      
system(system > 0) = psi_vecc((Ny+1):size(H,1)-size(cor_right,2));

C_coeff = psi_vecc(size(H,1)-size(cor_right,2)+1:end,1); % The last terms in psi_vec are the C coefficients! Same correction here, use proper psi_vec

psi_right = zeros(opening+extra_sites,res_length);

for y = 1:opening+extra_sites % Is there something wrong with this loop?!
    
    for x = (Nx+1):(Nx+res_length+1)
        
        sum = 0;
        
        for m = 1:opening+extra_sites
            
            sum = sum + C_coeff(m)*sin(pi*m*y*a/l_right)*exp(1i*k_prim(m)*x*a);
            
        end
        
        psi_right(y,(x-Nx-0)) = sum;
        
    end
    
end

psi = [zeros(extra_sites/2 + 1,Nx+2); system; zeros(extra_sites/2 + 1,Nx+2)];
psi_right = [zeros(1,res_length+1); psi_right; zeros(1,res_length+1)];
psi_full = [psi(:,1:(Nx + 1)),psi_right];

[d_x,d_y]=gradient(psi_full,a);

current_x = (h_bar/(2*m_e*1i))*(conj(psi_full).*d_x - psi_full.*conj(d_x));
current_y = (h_bar/(2*m_e*1i))*(conj(psi_full).*d_y - psi_full.*conj(d_y));

f11 = figure;
f22 = figure;

x = [0*1e9:a/1e-9:(Nx+1+res_length)*a/1e-9];
y = [-(opening+extra_sites+1)*(5/2)*a/5e-9:a/1e-9:(opening+extra_sites+1)*(5/2)*a/5e-9];

figure(f11);

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
line(0:x_lgth/10,y_lgth*ones(1,y_lgth/10+1),z_max*ones(1,y_lgth/10+1),'Color','red','LineWidth',2)
line(5:x_lgth,y_lgth:2*y_lgth-5,z_max*ones(1,y_lgth-4),'Color','red','LineWidth',2)
line(x_lgth*ones(1,6),2*y_lgth-5:2*y_lgth,z_max*ones(1,6),'Color','red','LineWidth',2)

line(0:x_lgth/10,-y_lgth*ones(1,y_lgth/10+1),z_max*ones(1,y_lgth/10+1),'Color','red','LineWidth',2)
line(5:x_lgth,-y_lgth:-1:-(2*y_lgth-5),z_max*ones(1,y_lgth-4),'Color','red','LineWidth',2)
line(x_lgth*ones(1,6),-(2*y_lgth-5):-1:-2*y_lgth,z_max*ones(1,6),'Color','red','LineWidth',2)

hold off

ax = gca;
ax.SortMethod;
ax.SortMethod = 'childorder';
% ax.SortMethod = 'depth';

xlabel('x (nm)')
ylabel('y (nm)')

zlabel('|\Psi|^2')
title(['Conic |\Psi|^2 for E = ' num2str(E*1e3) ' meV'])

ylim([-100 100])

figure(f22);

[X,Y] = meshgrid([0*1e9:a/1e-9:(Nx+1+res_length)*a/1e-9],[-(opening+extra_sites+1)*(5/2)*a/5e-9:a/1e-9:(opening+extra_sites+1)*(5/2)*a/5e-9]);

surf_plot_2 = surf(X,Y,sqrt(current_x.^2 + current_y.^2));

view(2)

hold on

% z_max = max(max(get(surf_plot_1,'Zdata')))
z_max = .1;
x_lgth = a*1e9*(Nx+1);
y_lgth = a*1e9*(Ny+1)/2;
line(0:x_lgth/10,y_lgth*ones(1,y_lgth/10+1),z_max*ones(1,y_lgth/10+1),'Color','red','LineWidth',2)
line(5:x_lgth,y_lgth:2*y_lgth-5,z_max*ones(1,y_lgth-4),'Color','red','LineWidth',2)
line(x_lgth*ones(1,6),2*y_lgth-5:2*y_lgth,z_max*ones(1,6),'Color','red','LineWidth',2)

line(0:x_lgth/10,-y_lgth*ones(1,y_lgth/10+1),z_max*ones(1,y_lgth/10+1),'Color','red','LineWidth',2)
line(5:x_lgth,-y_lgth:-1:-(2*y_lgth-5),z_max*ones(1,y_lgth-4),'Color','red','LineWidth',2)
line(x_lgth*ones(1,6),-(2*y_lgth-5):-1:-2*y_lgth,z_max*ones(1,6),'Color','red','LineWidth',2)

hold off

ax = gca;
ax.SortMethod;
ax.SortMethod = 'childorder';
% ax.SortMethod = 'depth';

xlabel('x (nm)')
ylabel('y (nm)')
zlabel('|J|')
title(['Conic |J| for E = ' num2str(E*1e3) ' meV'])

ylim([-100 100])





