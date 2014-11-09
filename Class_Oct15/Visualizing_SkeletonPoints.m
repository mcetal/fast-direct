% Plotting Ranks and Skeletonized Points
%-----------------------------------------

% Working on exterior multiply-connected ellipse problem

% This script applies the ID to the horizontal and vertical off diagonal
% blocks of the matrix M
% The k skeleton points chosen by the ID are plotted on the contours
%  compressed H^(i) : correspond to k rows, or skeleton source points that
%                     are chosen on gamma_i (ith section of contour)
%  compressed V^(i) : correspond to k columns, or skeleton target points
%                     that are chosen on gamma_i
% Replicating figures shown in Lecture 9 of Fast Direct Solvers Workshop
% 
% Natalia, Oct. 15,2014

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
title(['Exterior Domain' ' (N= ' num2str(N) ' points on each contour)']);
legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');

% ROW BLOCKS
H=[]; % Holds all blocks H^(i)
for i=1:4
    rowBlock=M((i-1)*N+1:(i-1)*N+N,:);
    rowBlock(1:N, (i-1)*N+1:(i-1)*N+N)=0;
    rowBlock(rowBlock==0)=[];
    rowBlock=reshape(rowBlock, [N,3*N]);
    H=[H;rowBlock];
end

% Plotting singular values of H^(i)
hFig=figure;
set(hFig, 'Position', [150 260  1000 400])
subplot(1,2,1);
colors={'b.','k.','r.','g.'};
for i=1:4
    s=svd(H((i-1)*N+1:(i-1)*N+N,:));
    semilogy(s,colors{i},'MarkerSize',5);
    hold on
end
title('Singular Values of H^{(i)}');

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
subplot(1,2,2);
colors={'b.','k.','r.','g.'};
for i=1:4
    s=svd(V(:,(i-1)*N+1:(i-1)*N+N));
    semilogy(s,colors{i},'MarkerSize',5);
    hold on 
end
    
legend('{\Gamma}_1','{\Gamma}_2','{\Gamma}_3','{\Gamma}_4','Location','Best');
title('Singular Values of V^{(i)}');

% SKELETON OF TARGETS (compression of V^(i))
acc=10^(-10);   % accuracy of compression
hFig=figure;
set(hFig, 'Position', [150 260  1000 400])
for i=1:4
    Vi=V(:, (i-1)*N+1:(i-1)*N+N);  % vertical block
    [T, I]=id_decomp(Vi,acc,'PGS');

    rank=size(T,1);
    Vics=Vi(:,I(1:rank));
    
    % For plotting off-diagonal blocks 
    % setting diagonal blocks to 0, vertical blocks to 1
    M2=zeros(size(M));
    M2(:,(i-1)*N+1:(i-1)*N+N)=1;  
    M2((i-1)*N+1:(i-1)*N+N,(i-1)*N+1:(i-1)*N+N)=0;
   
    subplot(1,2,1);
    
    % plotting vertical block
    spy(M2);
    title('Vertical Blocks V^{(i)}');
    
    % plotting k skeleton points
    subplot(1,2,2);
    for j=1:4
        if j==i
            plot(G{i}.rx(t(I(1:rank))),G{i}.ry(t(I(1:rank))),'b.','MarkerSize',5);
            hold on 
        else
            plot(G{j}.rx(t),G{j}.ry(t),'k.','MarkerSize',5);
            hold on 
        end
    end
    title(['Skeletonized target points with accuracy =' num2str(acc)]);
    xlabel('CLICK TO CONTINUE');
    w=waitforbuttonpress;
    if w==0
        disp('Button click');
    else 
        disp('Key press');
    end
    hold off; 
end

% SKELETON OF SOURCES (Compression of H^(i))
acc=10^(-10);
hFig=figure;
set(hFig, 'Position', [150 260  1000 400]);

for i=1:4
    Hi=H((i-1)*N+1:(i-1)*N+N,:);      % horizontal block
    [T, I]=id_decomp((Hi)',acc,'PGS');

    rank=size(T,1);
    Hics=Hi(I(1:rank),:);
    
    % For plotting off-diagonal blocks 
    % setting diagonal blocks to 0, horizontal blocks to 1
    M2=zeros(size(M));
    M2((i-1)*N+1:(i-1)*N+N,:)=1;  
    M2((i-1)*N+1:(i-1)*N+N,(i-1)*N+1:(i-1)*N+N)=0;

    subplot(1,2,1);
    
    % plotting horizontal block
    spy(M2,'r'); 
    title('Horizontal Blocks H^{(i)}');
    
    % plotting k skeleton points
    subplot(1,2,2);
    for j=1:4
        if j==i
            plot(G{i}.rx(t(I(1:rank))),G{i}.ry(t(I(1:rank))),'r.','MarkerSize',5);
            hold on 
        else
            plot(G{j}.rx(t),G{j}.ry(t),'k.','MarkerSize',5);
            hold on 
        end
    end
    title(['Skeletonized source points with accuracy =' num2str(acc)]);
    xlabel('CLICK TO CONTINUE');
    w=waitforbuttonpress;
    if w==0
        disp('Button click');
    else 
        disp('Key press');
    end
    hold off
end




