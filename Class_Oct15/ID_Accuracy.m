% Accuracy of Decomposition
%------------------------------
% Looking at accuracy of the ID on each off diagonal block, V^(i) and H^(i)
% Comparing each compressed block against the expected error bound given in
% Cheng et al. on pg. 1392

clear; clc; close all;
Nvec=[2.^(5:10)];  %Number of points on each ellipse N=32:1024
k=15; % FIXED RANK k

actual_err=zeros(length(Nvec),4);
exp_err=zeros(length(Nvec),4);
count=1;

% COL BLOCKS
for N=Nvec
    disp('Compressing Column Blocks');
    disp(['N=' num2str(N)]);
    [M rhs G t]=GetLinearSystemv2(N);
    
    V=[]; % Holds all blocks V^(i)
    
    for i=1:4
        colBlock=M(:,(i-1)*N+1:(i-1)*N+N);
        colBlock((i-1)*N+1:(i-1)*N+N, 1:N)=0;
        colBlock(colBlock==0)=[];
        colBlock=reshape(colBlock, [3*N,N]);
        V=[V colBlock];
    end
    
    for i=1:4
        Vi=V(:,(i-1)*N+1:(i-1)*N+N);
        
        [T,I]=id_decomp(Vi,k,'PGS');
        Ik=eye(k);
        Vcs=Vi(:,I(1:k));
        Decomp=Vcs*[Ik T];
        Id=eye(N);
        Id=Id(:,I);
        
        DecompVi=Decomp*Id';
        actual_err(count,i)=norm(Vi-DecompVi);
       
        s=svd(Vi);
        eps_sigma=s(k+1);
        exp_err(count,i)=eps_sigma*sqrt(1+k*(N-k));
        
    end
    count=count+1;
end

% Plotting actual and expected error of compression
colors={'b.','k.','r.','g.'};
colors2={'b-','k-','r-','g-'};
colors3={'b--','k--','r--','g--'};

for i=1:4
    semilogy(Nvec,actual_err(:,i),colors{i},'MarkerSize',5);
    hold on
    a=semilogy(Nvec,actual_err(:,i),colors2{i});
    b=semilogy(Nvec,exp_err(:,i),colors3{i});   
end
title(['Accuracy of ID on Vertical Blocks V^{(i)} Rank =' num2str(k)]);
legend([a b],'Actual Error', 'Expected Error','Location','Best');

% ROW BLOCKS
actual_err=zeros(length(Nvec),4);
exp_err=zeros(length(Nvec),4);
count=1;

for N=Nvec
    disp('Compressing Row Blocks');
    disp(['N=' num2str(N)]);
    [M rhs G t]=GetLinearSystemv2(N);
    H=[];     % Holds all blocks of H^(i)
    for i=1:4
        rowBlock=M((i-1)*N+1:(i-1)*N+N,:);
        rowBlock(1:N, (i-1)*N+1:(i-1)*N+N)=0;
        rowBlock(rowBlock==0)=[];
        rowBlock=reshape(rowBlock, [N,3*N]);
        H=[H;rowBlock];
    end
    
    for i=1:4
        Hi=(H((i-1)*N+1:(i-1)*N+N,:));
        [T,I]=id_decomp(Hi',k,'PGS');
        
        HiT=Hi';
        Hcs=HiT(:,I(1:k));
        
        IN=eye(N);
        IN=IN(:,I);
        Ik=eye(k);
        Decomp=[Ik;(T)']*Hcs';
        
        
        actual_err(count,i)=norm(Hi-IN*Decomp);
        
        s=svd(Hi);
        eps_sigma=s(k+1);
        exp_err(count,i)=eps_sigma*sqrt(1+k*(N-k));
        
    end
    count=count+1;
end

% Plotting actual and expected error of compression
figure;
colors={'b.','k.','r.','g.'};
colors2={'b-','k-','r-','g-'};
colors3={'b--','k--','r--','g--'};
for i=1:4
    
    semilogy(Nvec,actual_err(:,i),colors{i},'MarkerSize',5);
    hold on
    a=semilogy(Nvec,actual_err(:,i),colors2{i});
    b=semilogy(Nvec,exp_err(:,i),colors3{i});
end
title(['Accuracy of ID on Horizontal Blocks H^{(i)} Rank = ' num2str(k)]);
legend([a b],'Actual Error', 'Expected Error','Location','Best');
