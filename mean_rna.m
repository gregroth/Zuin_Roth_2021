%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mean_rna - This function calculate the mean number of mrna predicted by
%the enhancer-promoter model for a specific contact probability
%
%States
%           s1/off/close s2/off/close ... sn/off/close | s1/on/close ...
%           sn/on/close || s1/off/far ... sn/off/far | s1/on/far ...
%           sn/on/far
%
%Note  the mrna degradation rate delta = 1
%
%Input
%           - parm=[nbs,kfar,kfor,kback,konL,koffL,konH,koffH,muL,muH]
%           - cp= contact probability
%
%Output
%           - mean_rna =mean number of mrna
%
% Author: Gregory Roth
%
%   original version: 19.02.2020,
%   last version: 15.04.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_rna = mean_rna(parm,cp)
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
    for i=1:length(cp)
        kclose=kfar*cp(i)/(1-cp(i)); %kclose is calculated from kfar and cp
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
        
        
        
        
        Q=[BL11,BL12,BL13,BL14
            BL21,BL22,BL23,BL24
            BL31,BL32,BL33,BL34
            BL41,BL42,BL43,BL44];
        Q=-diag(sum(Q,2))+Q;
        
        %Calculate the stationary distribution (equation 8 in the SI)
        sdtemp=null(Q');
        sd=sdtemp./sum(sdtemp);%normalisation of the eigenvector
        [s1,s2]=size(sd);
         if s2>1
            sd=sd(:,end);
        end
        trates=[zeros(1,nbs),mudeltaL.*ones(1,nbs-1),mudeltaH,zeros(1,nbs),mudeltaL.*ones(1,nbs-1),mudeltaH]; %initiation rates, i.e. vector mu in SI
        mean_rna(i)=trates*sd; 
    end 
    mean_rna=mean_rna';
end

