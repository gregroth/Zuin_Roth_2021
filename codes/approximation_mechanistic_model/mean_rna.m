%% mean_rna -- this function calculate the mean number of mrna with the variable kon two-state model
%
%
%This function calculate the mean # rna as a function of cp and the
%paramters of the two state model
%
%
%the model
% 
%it is a two state model where the parameter kon depends on
%the contact probability cp through a hill function with constant c and exponent N
%
%Input
%           - totparam=vector with parameter values
%                       (N,konB,konE,koff,mu,beta)
%           - cp (contact probability)
%
%Output
%           - mean_rna = mean number of mrna 
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_rna = mean_rna(totparam,cp)
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

