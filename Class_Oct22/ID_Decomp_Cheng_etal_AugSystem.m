
% Compressing and Solving Augmented System using ID Decomposition from Cheng et al.
% Plots resulting matrix structure of non-augmented and augmented system

% Functions used: GetAugSystem.m, GetOffDiagBlock.m, id_decomp 

% Natalia Oct. 21, 2014

clear; clc; close all;
k=15; N=64;
[M rhs B C D]=GetAugSystem(N);
AugM=zeros(4*N+4);  
CompM=M;            % AugM and CompM updated as we go
R=[];               % holds block matrices R^(i)

for i=1:4
    % Get H^(i) using only M
    [Hi]=GetOffDiagBlock(CompM,'H',i,N);  
    [T,I]=id_decomp(Hi',15,'PGS');
    
    % building matrix L^(i)
    Ik=eye(k);
    Inmk=eye(N-k);
    IN=eye(N);
    IN=IN(:,I);
    Li=[Ik zeros(k,N-k); -T' Inmk]*IN';
    
    % compressed block with zeros introduced
    CompHi=Hi(I(1:k),:);
    CompHi=[CompHi;zeros(N-k,3*N)];
    
    % placing compressed block back into matrix and updating diagonal block
    CompM((i-1)*N+1:(i-1)*N+N,:)=[CompHi(:,1:(i-1)*N) Li*CompM((i-1)*N+1:(i-1)*N+N,(i-1)*N+1:(i-1)*N+N) CompHi(:,(i-1)*N+1:end)];
    AugM(1:4*N,1:4*N)=CompM;
    
    % updating L^(i)*Bi
    AugM((i-1)*N+1:(i-1)*N+N,4*N+1:end)=Li*B((i-1)*N+1:(i-1)*N+N,:);
    
    % updating  rhs: L^(i)*(-2f)
    rhs((i-1)*N+1:(i-1)*N+N)=Li*rhs((i-1)*N+1:(i-1)*N+N);
    
    % Get V^(i)
    [Vi]=GetOffDiagBlock(CompM,'V',i,N);
    [T,I]=id_decomp(Vi,15,'PGS');
    
    % building matrix R^(i)
    Ik=eye(k);
    IN=eye(N);
    Inmk=eye(N-k);
    IN=IN(:,I);
    Ri=IN*[Ik -T; zeros(N-k,k) Inmk];
    
    % compressed block with zeros introduced
    CompVi=Vi(:,I(1:k));
    CompVi=[CompVi zeros(3*N,N-k)];
    
    % placing compressed block back into matrix and updating diagonal block
    CompM(:,(i-1)*N+1:(i-1)*N+N)=[CompVi(1:(i-1)*N,:); CompM((i-1)*N+1:(i-1)*N+N,(i-1)*N+1:(i-1)*N+N)*Ri; CompVi((i-1)*N+1:end,:)];
    AugM(1:4*N,1:4*N)=CompM;
    
    % updating Ci*R^(i)
    AugM(4*N+1:end,(i-1)*N+1:(i-1)*N+N)=C(:,(i-1)*N+1:(i-1)*N+N)*Ri;
    
    % saving R^(i) for later
    R=[R;Ri];
end

% Adding on D to AugM
AugM(4*N+1:4*N+4,4*N+1:4*N+4)=D;

figure;
spy(CompM);
hold on
title('Compressed Matrix M, N=64');

figure;
spy(AugM);
hold on
title('Compressed Augmented Matrix AugM, N=64');

% Solve System
%-----------------
sigma=AugM\rhs;
R=blkdiag(R(1:N,:),R(N+1:2*N,:),R(2*N+1:3*N,:),R(3*N+1:4*N,:),eye(4));
sigma=R*sigma;

figure;
plot(sigma(1:4*N),'r');
title('Solution to Density Function from Compressed System, N=64');
