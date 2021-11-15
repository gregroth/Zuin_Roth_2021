%% MLtot_2state -- this function calculate likelihood of the FACS and FISH data with the variable kon two-state model
%
%
%Input
%           - freeparam=vector of 0 and 1 with 1 if the param is free
%           - fixedparam=vector with values of the fixed parameters
%           - x = vector with the values of the free parameter
%           - maxRna = vector with maximum # mrna per FISH clone
%           - dataVal = list of histogram value for each FISH clone
%           - binsize,trdata_l,cpdata_l,cpdata_fish
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


function out=MLtot_2state(freeparam,fixedparam,x,maxRna,dataVal,binsize,trdata_l,cpdata_l,cpdata_fish)
    totparam=zeros(1,6);
    totparam(freeparam==1)=x;
    totparam(freeparam==0)=fixedparam;
    N=length(trdata_l); %number of data points
    sigma=1; 

    part1=N*log(sigma*sqrt(2*pi))+sum((trdata_l-mean_rna_mixedparam(freeparam,fixedparam,x,cpdata_l)).^2)/(2*sigma^2);

    MLtemp=0;
    for i=1:length(maxRna)
        konL=totparam(1);
        konH=totparam(2);
        koff=totparam(3);
        mu=totparam(4);
        c=totparam(5);
        N=totparam(6);
        hillfactor=cpdata_fish(i)^N/(c+cpdata_fish(i)^N);
        kon=konL+hillfactor*(konH-konL);
        MLtemp=MLtemp+MLdist_2state([kon,koff,mu],1,maxRna(i),dataVal{i},binsize(i));
    end
    out=MLtemp+part1;
end
