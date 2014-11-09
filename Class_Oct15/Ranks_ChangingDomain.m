% Studying Ranks of Horizontal and Vertical Blocks for Changing Domains
%-----------------------------------------------------------------------
% Looking at effects of bringing contours closer together

% Natalia, Oct. 15

clear; clc; close all;
for k=1:3  % 3 different domains 
    N=64;
    [M rhs G t]=GetLinearSystemv3(N,k);
    
    % ROW BLOCKS
    H=[];     % Holds all blocks H^(i)
    for i=1:4
        rowBlock=M((i-1)*N+1:(i-1)*N+N,:);
        rowBlock(1:N, (i-1)*N+1:(i-1)*N+N)=0;
        rowBlock(rowBlock==0)=[];
        rowBlock=reshape(rowBlock, [N,3*N]);
        H=[H;rowBlock];
    end
    
    % Plotting singular values of H^(i)
    hFig=figure(2);
    set(hFig, 'Position', [100 20  600 600])
    subplot(3,2,2*(k-1)+1);
    colors={'b.','k.','r.','g.'};
    for i=1:4
        s=svd(H((i-1)*N+1:(i-1)*N+N,:));
        semilogy(s,colors{i},'MarkerSize',5);
        hold on
    end
    title(['Singular Values of H^{(i)}, Domain ' num2str(k)]);
    
    % COLUMN BLOCKS
    V=[];      % Holds all blocks V^(i)
    for i=1:4
        colBlock=M(:,(i-1)*N+1:(i-1)*N+N);
        colBlock((i-1)*N+1:(i-1)*N+N, 1:N)=0;
        colBlock(colBlock==0)=[];
        colBlock=reshape(colBlock, [3*N,N]);
        V=[V colBlock];
    end
    
    % Plotting singular values of V^(i)
    subplot(3,2,2*k);
    colors={'b.','k.','r.','g.'};
    for i=1:4
        s=svd(V(:,(i-1)*N+1:(i-1)*N+N));
        semilogy(s,colors{i},'MarkerSize',5);
        hold on
    end
    
   % legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');
    title(['Singular Values of V^{(i)}, Domain ' num2str(k)]);
end
