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
%           - cp (contact probability)
%           - parm=[kon,koff,mu,c,N], cp
%
%Output
%           - meanrnaout = mean number of mrna 
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meanrnaout = mean_rna(parm,cp)
konL=parm(1);
konH=parm(2);
koff=parm(3);
mu=parm(4);
c=parm(5);
N=parm(6);
    for i=1:length(cp)
        hillfactor=cp(i)^N/(c+cp(i)^N);
        kontemp=konL+hillfactor*(konH-konL);
        kofftemp=koff;
        mutemp=mu;
        meanrna(i)=kontemp/(kontemp+kofftemp)*mutemp;
    end 
    meanrnaout=meanrna';
end

