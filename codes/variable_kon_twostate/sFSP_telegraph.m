%% sFSP_telegraph -- this function calculates an approximation of the steady state probability distribution of the two-state model. The algorithm follows the stationary Finite State Projectio
%
%Note: the approximation method is provided by Ankit Gupta, Jan Mikelson and Mustafa in
%"A ﬁnite state projection algorithm for the stationary solution of the
%chemical master equation" 
%
%Input
%
%           - parm=[kon,koff,kini,delta], parameter vector
%           - maxRna = maximum number of RNA molecules
%Output
%
%           -  dist.rna= vector containing the marginal distribution of the
%              number of RNA
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function distRNA = sFSP_telegraph(parm,maxRna) 
    kon=parm(1);
    koff=parm(2);
    mu=parm(3);
    delta=parm(4);
    
    S=[-koff,koff
    kon,-kon];

    dim=[2,maxRna]; %we expend the size of the last dimension to obtain at the end the correct diagonal
   
    ee=ones(1,maxRna-1);
    el=1:1:maxRna-1;
    trstates=[mu,0];
    Ttemp=diag(trstates);
   
    Qtemp1=kron(eye(maxRna),S)+kron(diag(ee,1),Ttemp)+kron(diag(el,-1),delta.*eye(2)); 
    
    cn=zeros(prod(dim),1); %vector having the sum of the rates of the outgoing transition
    cn(end-1:end)=trstates';
    bl=zeros(1,length(cn));
    bl(end-1)=1;
    Qabs=cn*bl;
    
    Qtemp2=Qtemp1+Qabs;
    Q=sparse(Qtemp2-diag(sum(Qtemp2,2)));
    
    [V,D] = eigs(Q',1,10^-2);
    disttot=abs(V./sum(V));
    %calculate the marginal distribution (i.e. dist. of # RNA)
    mattemp=kron(eye(maxRna),ones(1,2)); 
    distRNA=mattemp*disttot;
end