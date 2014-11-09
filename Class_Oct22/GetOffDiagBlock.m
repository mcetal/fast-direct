function [Block]=GetOffDiagBlock(M, dir,index,N)
%------------------------------------------------------
% INPUTS: M : 4*N x 4*N system
%        dir: string specifying 'H' or 'V'
%      index: index of block 
% OUTPUTS: horizontal or vertical off-diagonal block 
%------------------------------------------------------

if strcmp(dir,'H') % Get Horizontal Block
    Block=M((index-1)*N+1:(index-1)*N+N,:);
    Block(1:N,(index-1)*N+1:(index-1)*N+N)=NaN;
    Block(isnan(Block))=[];
    Block=reshape(Block,[N, 3*N]);    
else  % Get Vertical Block
    Block=M(:,(index-1)*N+1:(index-1)*N+N);
    Block((index-1)*N+1:(index-1)*N+N,1:N)=NaN;
    Block(isnan(Block))=[];
    Block=reshape(Block, [3*N,N]);
end 
