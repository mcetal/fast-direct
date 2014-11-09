
% Exterior Dirichlet Problem
%-----------------------------

% This script solves Laplace's equation for the exterior Dirichlet problem. 
% Domain : unbounded region exterior to 4 well-separated ellipses.
% Using the method outlined in Greengard, Greenbaum, and McFadden (1993).
%
% Natalia, Oct.14, 2014

% Defining 4 Ellipses
%----------------------
clear; clc; close all;

a1=0.5; b1=0.75;   % major and minor axes
a2=0.4; b2=0.25;
a3=0.25; b3=0.5;
a4=0.5;  b4=0.25;

cx1=0; cy1=0;       % centers
cx2=3.5; cy2=0.4;
cx3=0.2; cy3=-5.25;
cx4=4.5; cy4=-5.25;

per=10^(-15);                % perturbing slightly so centers don't correspond with grid
sx1=cx1+per; sy1=cy1+per;       
sx2=cx2+per; sy2=cy2+per;
sx3=cx3+per; sy3=cy3+per;
sx4=cx4+per; sy4=cy4+per;

sx=[sx1 sx2 sx3 sx4];
sy=[sy1 sy2 sy3 sy4];

% A structure is used to store information about each curve
% i.e. parametric represenation, normal vectors, etc. 
G{1}=MakeEllipse(a1,b1,cx1,cy1);  
G{2}=MakeEllipse(a2,b2,cx2,cy2);
G{3}=MakeEllipse(a3,b3,cx3,cy3);
G{4}=MakeEllipse(a4,b4,cx4,cy4);

% Number of points on each boundary
N=64;
h=2*pi/N;
t=(0:N-1)*h;

%Plotting Domain: 
figure;
colors={'b.','k.','r.','g.'};
for i=1:4
    plot(G{i}.rx(t),G{i}.ry(t),colors{i},'MarkerSize',5);
    hold on 
end
title('Exterior Domain');
legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');

% System for Integral Equation:
%--------------------------------

% combining information altogether for all ellipses
rx=[]; ry=[]; nx=[]; ny=[];
for i=1:4
    rx=[rx G{i}.rx(t)];
    ry=[ry G{i}.ry(t)];
    nx=[nx G{i}.nx(t)];
    ny=[ny G{i}.ny(t)];
end

length=[]; curv=[];
for i=1:4
    for j=1:N
        Gilength(j)=G{i}.length(t(j));
        Gicurv(j)=G{i}.curv(t(j));
    end
    length=[length Gilength];
    curv=[curv Gicurv];
end

% Exact Solution (and Boundary Data)
%f=@(x,y) exp(x./(x.^2+(y).^2)).*cos((y)./(x.^2+(y).^2));
f=@(x,y) real(1./((x+1i*y)-(sx(1)+1i*sy(1)))+1./((x+1i*y)-(sx(2)+1i*sy(2)))+...
    1./((x+1i*y)-(sx(3)+1i*sy(3)))+1./((x+1i*y)-(sx(4)+1i*sy(4))));

% formula for kernel (grad log * n)
kernel=@(nx,ny,rx1,rx2,ry1,ry2,length) (nx*(rx2-rx1)+ny*(ry2-ry1))./(((rx1-rx2).^2+(ry1-ry2).^2)*length);

K=zeros(4*N,4*N); % matrix representing kernel
for i=1:4*N
    for j=1:4*N
        if i==j  
            K(i,j)=h*((1/2)*curv(i)*length(i)+1);
        else 
            K(i,j)=h*(kernel(nx(j),ny(j),rx(i),rx(j),ry(i),ry(j),length(j))*length(j)+1);
        end
    end
end
    
B=zeros(4*N,4);
for i=1:4*N
    for j=1:4
        B(i,j)=-2*log(sqrt((rx(i)-sx(j)).^2+(ry(i)-sy(j)).^2));
    end
end

I=eye(4*N,4*N);
rhs=[-2*f(rx, ry)'; zeros(4,1)];

%solving for sigma
M=(I-(1/pi)*K);

hvec=h*ones(1,N);
C=blkdiag(hvec,hvec,hvec,zeros(1,N));

D=[zeros(3,4); ones(1,4)];

M=[I-(1/pi)*K B; C D];
sigma=M\rhs;          % Solution for sigma and coefficients A
A=sigma(end-3:end);
sigma=sigma(1:end-4);

% Reconstruct Solution
%-----------------------
g1=-6:0.0123:6;        % grid for evaluating solution
g2=-6:0.0123:6;

[x y] = meshgrid(g1,g1);     
        
u=zeros(size(x));
for i=1:4*N
    u=u+(1/(2*pi))*h*(kernel(nx(i),ny(i),x,rx(i),y,ry(i),length(i))*length(i)+1)*sigma(i);
    for j=1:4
       u=u-2*A(j)*log(sqrt((x-sx(j)).^2+(y-sy(j)).^2));
    end
        
end

figure;   % Error Plot
imagesc(g1,g2,log10(abs(u-f(x,y)))); caxis([-16 0]); colorbar;
title(['Exterior Dirichlet Laplace log error, N=' num2str(N)] );
axis xy equal tight;
hold on
plot(rx,ry,'k.','markersize',5);
axis xy equal tight;

% Convergence Rate of Error
%----------------------------
Nvec=5:5:110;
count=1;
err=zeros(size(Nvec));

for N=Nvec
    disp(['N=' num2str(N) ' out of ' num2str(Nvec(end))]);
    
    h=2*pi/N;
    t=(0:N-1)*h;
    
    % System for Integral Equation:
    %--------------------------------
    % combining information altogether for all ellipses
    rx=[]; ry=[]; nx=[]; ny=[];
    for i=1:4
        rx=[rx G{i}.rx(t)];
        ry=[ry G{i}.ry(t)];
        nx=[nx G{i}.nx(t)];
        ny=[ny G{i}.ny(t)];
    end

    length=[]; curv=[];
    for i=1:4
        for j=1:N
            Gilength(j)=G{i}.length(t(j));
            Gicurv(j)=G{i}.curv(t(j));
        end
        length=[length Gilength];
        curv=[curv Gicurv];
    end

    % Exact Solution (and Boundary Data)
    %f=@(x,y) exp(x./(x.^2+(y).^2)).*cos((y)./(x.^2+(y).^2));
    f=@(x,y) real(1./((x+1i*y)-(sx(1)+1i*sy(1)))+1./((x+1i*y)-(sx(2)+1i*sy(2)))+...
        1./((x+1i*y)-(sx(3)+1i*sy(3)))+1./((x+1i*y)-(sx(4)+1i*sy(4))));

    % formula for kernel (grad log * n)
    kernel=@(nx,ny,rx1,rx2,ry1,ry2,length) (nx*(rx2-rx1)+ny*(ry2-ry1))./(((rx1-rx2).^2+(ry1-ry2).^2)*length);

    K=zeros(4*N,4*N); % matrix representing kernel
    for i=1:4*N
        for j=1:4*N
            if i==j  
                K(i,j)=h*((1/2)*curv(i)*length(i)+1);
            else 
                K(i,j)=h*(kernel(nx(j),ny(j),rx(i),rx(j),ry(i),ry(j),length(j))*length(j)+1);
            end
        end
    end
    
    B=zeros(4*N,4);
    for i=1:4*N
        for j=1:4
            B(i,j)=-2*log(sqrt((rx(i)-sx(j)).^2+(ry(i)-sy(j)).^2));
        end
    end

    I=eye(4*N,4*N);
    rhs=[-2*f(rx, ry)'; zeros(4,1)];

    M=(I-(1/pi)*K);

    hvec=h*ones(1,N);
    C=blkdiag(hvec,hvec,hvec,zeros(1,N));

    D=[zeros(3,4); ones(1,4)];

    M=[I-(1/pi)*K B; C D];
    sigma=M\rhs;           % Solution for sigma and coefficients A
    A=sigma(end-3:end);
    sigma=sigma(1:end-4);

    % Reconstruct Solution (Well away from boundary)
    
    g1=-10:0.0123:-9;
    g2=-10:0.0123:-9;

    [x y] = meshgrid(g1,g1);     
        
    u=zeros(size(x));
    for i=1:4*N
        u=u+(1/(2*pi))*h*(kernel(nx(i),ny(i),x,rx(i),y,ry(i),length(i))*length(i)+1)*sigma(i);
        for j=1:4
            u=u-2*A(j)*log(sqrt((x-sx(j)).^2+(y-sy(j)).^2));
        end
    end
    
    err(count)=norm(u-f(x,y),inf);
    count=count+1;
end

figure; % Plotting error
semilogy(Nvec,err,'.','MarkerSize',5);
hold on 
semilogy(Nvec,err);
title('Convergence Rate of Error (Well Away from Boundary)');
ylabel('L_{\infty} error');
xlabel('N');

