%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fsp_ep - This function calculates an approximation of the steady-state 
%probability distribution of the enhancer-promoter model
%
% 
%The algorithm follows the Finite State Projectio
%approximation method provided by Munsky and Khammash, 2006
%
%States
%           s1/off/close s2/off/close ... sn/off/close | s1/on/close ...
%           sn/on/close || s1/off/far ... sn/off/far | s1/on/far ...
%           sn/on/far
%
%Note  the mrna degradation rate delta = 1
%
%Input
%
%           - parm=[nbs,kfar,kfor,kback,konL,koffL,konH,koffH,muL,muH], parameter vector
%           - maxRna = maximum number of RNA molecules
%           - nbs, number of waiting states
%           - p0 = initial distribution, vector of size maxRna+1
%           - t= time at which we calculate the approximation
%
%Output
%
%           -  dist.join= vector containing the distribution of the state
%              of the promoter, the state of the enhancer, and the number of RNA
%           -  dist.rna= vector containing the marginal distribution of the
%              number of RNA
%           -  dist.mean= mean of the # RNA
%           -  dist.var= variance of the # RNA
%
% Other m-files required: rate_matrix.m
%
% Author: Gregory Roth
%
%   original version: 19.02.2020,
%   last version: 15.04.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist = fsp_ep(parm,cp,maxRna,t)
if cp==1
cp=.999;
end
    nbs=ceil(parm(1));%number of regulatory steps
    p0=zeros(4*nbs*maxRna,1);%initial distribution
    p0(1)=1;
    QQ=sparse(rate_matrix(maxRna,parm,cp));%calculate the truncated rate matrix
   
    disttemp=expm(QQ'*t)*p0;%calculate the distribution at time t
    %calculate the marginal distribution (i.e. dist. of # RNA)
    mattemp=kron(eye(maxRna),ones(1,4*nbs)); 
    distRNA=mattemp*disttemp;
    dist.join=disttemp;
    dist.rna=distRNA;
    
    %mean # of RNA
    ee=(0:1:maxRna-1);
    dist.mean=ee*distRNA;
    
    %variance of # of RNA
    dist.var=ee.^2*distRNA -(ee*distRNA).^2;
end