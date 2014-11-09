%--------------------------------------------------------------------------
% This script compares the solutions to the original augmented system and 
% compressed system.
% Solutions are compared for fixed N across varying k, and for fixed k, 
% across varying N
% The original and cormpressed solutions come from the functions: 
%  GetSolution_OrigSystem(N)
%  GetSolution_CompSystem(N,k)

% Natalia, Oct.21, 2014
%--------------------------------------------------------------------------

clear; clc; close all;

% For Fixed N, Varying k
N=64; 
kvec=[15:64];
diff=zeros(size(kvec,1),1);
count=1;

disp('Fixed N, Varying k');
for k=kvec
    disp(['k=' num2str(k) ' out of ' num2str(kvec(end))]);
    sigma_orig=GetSolution_OrigSystem(N);
    sigma_comp=GetSolution_CompSystem(N,k);
    diff(count)=norm(sigma_orig-sigma_comp,inf);
    count=count+1;
end

figure;
semilogy(kvec,diff);
hold on 
semilogy(kvec,diff,'.','MarkerSize',5);
title('Difference between Original and Compressed Density for varying k, N=64');
xlabel('k');
ylabel('L_{\infty} norm');

figure;
plot(sigma_comp);
hold on
plot(sigma_orig,'r');
title('Approximate density \sigma');
legend('Compressed', 'Original','Location','Best');

% For Fixed k, Varying N
Nvec=[16:5:100];
k=15;
diff=zeros(size(kvec,1),1);
count=1;

disp('Fixed k, Varying N');
for N=Nvec
    disp(['N=' num2str(N) ' out of ' num2str(Nvec(end))]);
    sigma_orig=GetSolution_OrigSystem(N);
    sigma_comp=GetSolution_CompSystem(N,k);
    diff(count)=norm(sigma_orig-sigma_comp,inf);
    count=count+1;
end

figure;
semilogy(Nvec,diff);
hold on 
semilogy(Nvec,diff,'.','MarkerSize',5);
title('Difference between Original and Compressed Density for varying N, k=15');
xlabel('k');
ylabel('L_{\infty} norm');

% RECONSTRUCTING SOLUTION in exterior domain
%---------------------------------------------
% Reconstructing solution using sigma from compressed system to compare
% against solution u

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

Nvec=64:5:130;
count=1;
kvec=[15 20 64];
err_orig=zeros(size(Nvec));
err_comp1=zeros(size(Nvec));
err_comp2=zeros(size(Nvec));
err_comp3=zeros(size(Nvec));

disp('Reconstructing Solution from density of compressed and uncompressed system');
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
    
    sigma_orig=GetSolution_OrigSystem(N);
    countk=1;
    sigma_comp=[];
    for k=kvec
        sigma_comp=[sigma_comp GetSolution_CompSystem(N,k)];
        countk=countk+1;
    end
    
    A_orig=sigma_orig(end-3:end);
    for j=1:size(kvec,2)
        A_comp(:,j)=sigma_comp(end-3:end);
    end
    sigma_comp(end-3:end,:)=[];
    sigma_orig=sigma_orig(1:end-4);
   
    % Reconstruct Solution (Well away from boundary)
    
    g1=-10:0.0123:-9;
    g2=-10:0.0123:-9;

    [x y] = meshgrid(g1,g1);     
        
    u_orig=zeros(size(x));
    u_comp1=zeros(size(x));
    u_comp2=zeros(size(x));
    u_comp3=zeros(size(x));
    for i=1:4*N
        u_orig=u_orig+(1/(2*pi))*h*(kernel(nx(i),ny(i),x,rx(i),y,ry(i),length(i))*length(i)+1)*sigma_orig(i);
        u_comp1=u_comp1+(1/(2*pi))*h*(kernel(nx(i),ny(i),x,rx(i),y,ry(i),length(i))*length(i)+1)*sigma_comp(i,1);
        u_comp2=u_comp2+(1/(2*pi))*h*(kernel(nx(i),ny(i),x,rx(i),y,ry(i),length(i))*length(i)+1)*sigma_comp(i,2);
        u_comp3=u_comp3+(1/(2*pi))*h*(kernel(nx(i),ny(i),x,rx(i),y,ry(i),length(i))*length(i)+1)*sigma_comp(i,3);
        for j=1:4
            u_orig=u_orig-2*A_orig(j)*log(sqrt((x-sx(j)).^2+(y-sy(j)).^2));
            u_comp1=u_comp1-2*A_comp(j,1)*log(sqrt((x-sx(j)).^2+(y-sy(j)).^2));
            u_comp2=u_comp2-2*A_comp(j,2)*log(sqrt((x-sx(j)).^2+(y-sy(j)).^2));
            u_comp3=u_comp3-2*A_comp(j,3)*log(sqrt((x-sx(j)).^2+(y-sy(j)).^2));
        end
    end
    
    err_orig(count)=norm(u_orig-f(x,y),inf);
    err_comp1(count)=norm(u_comp1-f(x,y),inf);
    err_comp2(count)=norm(u_comp2-f(x,y),inf);
    err_comp3(count)=norm(u_comp3-f(x,y),inf);
    count=count+1;
end

figure; % Plotting error
a=semilogy(Nvec,err_orig,'.','MarkerSize',5);
hold on 
semilogy(Nvec,err_orig);
title('Error in Solution (Well Away from Boundary)');
ylabel('L_{\infty} error');
xlabel('N');
hold on 

b=semilogy(Nvec,err_comp1,'r.','MarkerSize',5);
hold on 
semilogy(Nvec,err_comp1,'r');
title('Error in Solution');
ylabel('L_{\infty} error');
xlabel('N');

c=semilogy(Nvec,err_comp2,'k.','MarkerSize',5);
hold on 
semilogy(Nvec,err_comp2,'k');
title('Error in Solution');
ylabel('L_{\infty} error');
xlabel('N');

d=semilogy(Nvec,err_comp3,'g.','MarkerSize',5);
hold on 
semilogy(Nvec,err_comp3,'g');
title('Error in Solution ');
ylabel('L_{\infty} error');
xlabel('N');

legend([a,b,c,d], 'orig', ['k=' num2str(kvec(1))], ['k=' num2str(kvec(2))], ['k=' num2str(kvec(3))]);








