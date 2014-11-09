function [M rhs G t]=GetLinearSystemv2(N)

% This function produces a linear system for the solution of Laplace's
% equation for the exterior Dirichlet problem. 
% Domain : unbounded region exterior to 4 well-separated ellipses.
% Outputs: M: matrix containing the interactions between the boundaries of
%             each ellipse
%          rhs: right hand side vector containing boundary data
% If the system is solved we obtain the density function sigma: sigma=M\f
% Boundary Data is taken from the exact solution: 
%     u(x,y)=exp(x/(x^2+y^2)*cos(y/(x^2+y^2)
%
% Natalia, Oct.5, 2014

% Defining 4 Ellipses
%----------------------

a1=0.5; b1=0.75;   % major and minor axes
a2=0.4; b2=0.25;
a3=0.25; b3=0.5;
a4=0.5;  b4=0.25;

cx1=0; cy1=0;       % centers
cx2=3.5; cy2=0.4;
cx3=0.2; cy3=-5.25;
cx4=4.5; cy4=-5.25;

% A structure is used to store information about each curve
% i.e. parametric represenation, normal vectors, etc. 
G{1}=MakeEllipse(a1,b1,cx1,cy1);  
G{2}=MakeEllipse(a2,b2,cx2,cy2);
G{3}=MakeEllipse(a3,b3,cx3,cy3);
G{4}=MakeEllipse(a4,b4,cx4,cy4);

% Number of points on each boundary
h=2*pi/N;
t=(0:N-1)*h;

% Plotting Domain: 
% figure;
% colors={'b.','k.','r.','g.'};
% for i=1:4
%     plot(G{i}.rx(t),G{i}.ry(t),colors{i},'MarkerSize',5);
%     hold on 
% end
% title('Exterior Domain');
% legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');

% System for Integral Equation:
%--------------------------------

% combining information altogether for all ellipses
rx=[];
ry=[];
nx=[];
ny=[];
for i=1:4
    rx=[rx G{i}.rx(t)];
    ry=[ry G{i}.ry(t)];
    nx=[nx G{i}.nx(t)];
    ny=[ny G{i}.ny(t)];
end

length=[];
curv=[];
for i=1:4
    
    for j=1:N
        Gilength(j)=G{i}.length(t(j));
        Gicurv(j)=G{i}.curv(t(j));
    end
    length=[length Gilength];
    curv=[curv Gicurv];
end

% Exact Solution (and Boundary Data)
f=@(x,y) exp(x./(x.^2+(y).^2)).*cos((y)./(x.^2+(y).^2));

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
    
I=eye(4*N,4*N);
rhs=2*f(rx, ry)';

M=(-I+(1/pi)*K);


