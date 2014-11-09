% Plotting ranks of off diagonal blocks for varying values of N
%----------------------------------------------------------------

% Working on exterior multiply-connected ellipse problem

% Should confirm that the rank stays the same as N increases
% Natalia, Oct. 15

clear; clc; close all;
N=64;
[M rhs G t]=GetLinearSystemv2(N);

% plotting domain
figure;
colors={'b.','k.','r.','g.'};
for i=1:4
    plot(G{i}.rx(t),G{i}.ry(t),colors{i},'MarkerSize',5);
    hold on
end
title('Exterior Domain');
legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');

% ROW BLOCKS
Nvec=[64,128,256];
count=1;
figure;
for N=Nvec
    [M rhs G t]=GetLinearSystemv2(N);
    H=[]; % Holds all blocks H^(i)
    for i=1:4
        rowBlock=M((i-1)*N+1:(i-1)*N+N,:);
        rowBlock(1:N, (i-1)*N+1:(i-1)*N+N)=0;
        rowBlock(rowBlock==0)=[];
        rowBlock=reshape(rowBlock, [N,3*N]);
        H=[H;rowBlock];
    end
    
    % Plotting singular values of H^(i)
    subplot(3,2,(count-1)*2+1);
    colors={'b.','k.','r.','g.'};
    for i=1:4
        s=svd(H((i-1)*N+1:(i-1)*N+N,:));
        semilogy(s,colors{i},'MarkerSize',5);
        hold on
    end
    title(['Ranks of H^{(i)} N =' num2str(N)]);
    
    % COLUMN BLOCKS
    V=[]; % Holds all blocks V^(i)
    
    for i=1:4
        colBlock=M(:,(i-1)*N+1:(i-1)*N+N);
        colBlock((i-1)*N+1:(i-1)*N+N, 1:N)=0;
        colBlock(colBlock==0)=[];
        colBlock=reshape(colBlock, [3*N,N]);
        V=[V colBlock];
    end
    
    % Plotting singular values of V^(i)
    subplot(3,2,2*count);
    colors={'b.','k.','r.','g.'};
    for i=1:4
        s=svd(V(:,(i-1)*N+1:(i-1)*N+N));
        semilogy(s,colors{i},'MarkerSize',5);
        hold on
    end
    
    %legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');
    title(['Ranks of V^{(i)} N=' num2str(N)] );
    count=count+1;
end