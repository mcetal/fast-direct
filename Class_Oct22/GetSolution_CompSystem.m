function [sigma]=GetSolution_CompSystem(N,k)

% This function compresses the augmented system given from the function
% GetAugSystem using the ID Decomposition.

% The algorithm is based on the one outlined in Cheng et al. where the
% system is being left and right multiplied by matrices L^(i) and R^(i)
% The compressed system is then solved for density function sigma

% Inputs: N - number of points on each ellipse
%         k - rank of off-diagonal compression
% Outputs: sigma - approximate solution to density function and
%                  coefficients A1,..,A4
% Functions used: GetAugSystem, GetOffDiagBlock, id_decomp

[M rhs B C D]=GetAugSystem(N);
AugM=zeros(4*N+4);  
CompM=M;            % AugM and CompM updated as we go
R=[];               % holds matrices R^(i)

for i=1:4
    % Get H^(i) using only M
    [Hi]=GetOffDiagBlock(CompM,'H',i,N);  
    [T,I]=id_decomp(Hi',k,'PGS');
    
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
    
    % updating rhs: L^(i)*(-2f)
    rhs((i-1)*N+1:(i-1)*N+N)=Li*rhs((i-1)*N+1:(i-1)*N+N);
    
    % Get V^(i)
    [Vi]=GetOffDiagBlock(CompM,'V',i,N);
    [T,I]=id_decomp(Vi,k,'PGS');
    
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

% SOLVING SYSTEM
sigma=AugM\rhs;
R=blkdiag(R(1:N,:),R(N+1:2*N,:),R(2*N+1:3*N,:),R(3*N+1:4*N,:),eye(4));
sigma=R*sigma;
