%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cv_rna - This function builds the truncated rate matrix for the 
%enhancer-promoter model for a specific contact probability
%
%States
%           s1/off/close s2/off/close ... sn/off/close | s1/on/close ...
%           sn/on/close || s1/off/far ... sn/off/far | s1/on/far ...
%           sn/on/far
%
%Note  the mrna degradation rate delta = 1
%
%Input
%           - maxRNA = upper bound for the number of mRNA
%           - parm=[nbs,kfar,kfor,kback,konL,koffL,konH,koffH,muL,muH]
%           - cp= contact probability
%
%Output
%           - QQ =truncated rate matrix
%
% Author: Gregory Roth
%
%   original version: 19.02.2020,
%   last version: 15.04.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function QQ = rate_matrix(maxRna,parm,cp)
nbs=ceil(parm(1));
kfar=parm(2);
kfor=parm(3);
kback=parm(4);
konL=parm(5);
koffL=parm(6);
konH=parm(7);
koffH=parm(8);
mudeltaL=parm(9);
mudeltaH=parm(10);
delta=1;
        kclose=kfar*cp/(1-cp);%kclose is calculated from kfar and cp
        %calculation of the hyper-promoter-state transition matrix (i.e.
        %matrix K in the Supplementary Information)
        BL11=diag(kfor*ones(nbs-1,1),1)+diag(kback*ones(nbs-1,1),-1);
        BL12=konL.*eye(nbs);
        BL12(end,end)=konH;
        BL13=kfar.*eye(nbs);
        BL14=zeros(nbs);
        
        BL21=koffL*eye(nbs);
        BL21(end,end)=koffH;
        BL22=diag(kfor*ones(nbs-1,1),1)+diag(kback*ones(nbs-1,1),-1);
        BL23=zeros(nbs);
        BL24=kfar.*eye(nbs);
        
        BL31=kclose.*eye(nbs);
        BL32=zeros(nbs);
        BL33=diag(kback*ones(nbs-1,1),-1);
        BL34=BL12;
        
        BL41=zeros(nbs);
        BL42=kclose.*eye(nbs);
        BL43=BL21;
        BL44=diag(kback*ones(nbs-1,1),-1);
        
        
        
        
        S=[BL11,BL12,BL13,BL14
            BL21,BL22,BL23,BL24
            BL31,BL32,BL33,BL34
            BL41,BL42,BL43,BL44];
        S=-diag(sum(S,2))+S;
      
    dim=[4*nbs,maxRna+1]; %we expend the size of the last dimension to obtain at the end the correct diagonal
   
    ee=ones(1,maxRna);
    el=1:1:maxRna;
    trstates=[zeros(1,nbs),mudeltaL.*ones(1,nbs-1),mudeltaH,zeros(1,nbs),mudeltaL.*ones(1,nbs-1),mudeltaH];
    Ttemp=diag(trstates);
    
    %combine the hyper promoter state transition matrix witht the mRNA
    %production and degradation rate matrices
    QQtemp1=kron(eye(maxRna+1),S)+kron(diag(ee,1),Ttemp)+kron(diag(el,-1),delta.*eye(4*nbs)); 

    QQtemp2=diag(sum(QQtemp1,2));%matrix with the rate sum for each row on the diagonal
    QQtemp2=QQtemp1-QQtemp2;
    
    dim(end)=dim(end)-1; % remove the expension in the last dimension
    QQ=QQtemp2(1:prod(dim),1:prod(dim));
end

