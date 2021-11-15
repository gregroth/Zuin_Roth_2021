%% mean_rna_mixedparam -- this function calculate the mean number of mrna with the apparent two-state model
%
%
%This function calculate the mean # rna as a function of cp and the
%paramters of the apparent two state model
%
%
%the model
% 
%it is a two state model where the parameter kon depends on
%the contact probability cp through the nonlinear function described in eq. (31) in the suuplementary model desription
%
%Input
%           - freeparam=vector of 0 and 1 with 1 if the param is free
%           - fixedparam=vector of the fixed parameters
%           - cp (contact probability)
%           - param=value of the free parameters (among [N,konB,konE,koff,mu,beta])
%
%Output
%           - mean_rna = mean number of mrna 
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_rna = mean_rna_mixedparam(freeparam,fixedparam,x,cp)
    totparam=zeros(1,6);
    totparam(freeparam==1)=x;
    totparam(freeparam==0)=fixedparam;
    N=totparam(1);
    konB=totparam(2);
    konE=totparam(3);
    koff=totparam(4);
    mu=totparam(5);
    beta=totparam(6);
    for i=1:length(cp)
        alpha=cp(i).*beta;
        Hfactor=alpha^N*(1-alpha)/(1-alpha^(N+1));
        kon_app=konB+Hfactor*(konE-konB);
        mean_rna(i)=kon_app/(kon_app+koff)*mu;
    end 
    mean_rna=mean_rna';
end

