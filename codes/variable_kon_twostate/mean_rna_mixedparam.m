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
%the contact probability cp through a hill function with constant c and
%exponent N (see equation (3) in the supplementary model description).
%
%Input
%           - freeparam=vector of 0 and 1 with 1 if the param is free
%           - fixedparam=vector of the fixed parameters
%           - cp (contact probability)
%           - param=value of the free parameters (among [kon0,kon1,koff,mu,c,N])
%
%Output
%           - mean_rna = mean number of mrna 
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mean_rna = mean_rna_mixedparam(freeparam,fixedparam,param,cp)
totparam=zeros(1,6);
totparam(freeparam==1)=param;
totparam(freeparam==0)=fixedparam;
kon0=totparam(1);
kon1=totparam(2);
koff=totparam(3);
mu=totparam(4);
c=totparam(5);
N=totparam(6);
    for i=1:length(cp)
        hillfactor=cp(i)^N/(c+cp(i)^N);
        kontemp=kon0+hillfactor*(kon1-kon0);
        kofftemp=koff;
        mutemp=mu;
        mean_rna(i)=kontemp/(kontemp+kofftemp)*mutemp;
    end 
    mean_rna=mean_rna';
end

