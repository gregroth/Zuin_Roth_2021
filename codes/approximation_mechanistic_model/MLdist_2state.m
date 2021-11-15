%% MLdist_2state -- this function calculate likelihood of the FISH data with the variable kon two-state model
%
%
%Input
%           - x = vector of paramter (kon, koff, kini)
%           - delta = paramter delta (degradation rate)
%           - maxRna = cut-off of the number of RNA
%           - dataHist = histogram of the data (vector)
%           - Bwidth =  width of the bin
%
%Output
%           - like = maximum likelihood
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function like = MLdist_2state(x,delta,maxRna,dataHist,Bwidth)  
    parm=[x,delta];
    %ATTENTION this is to remove extremely small probabilities and avoid
    %ML=inf
    %remove negative probabilities
    dataHist(dataHist<0)=0;
    distRNA=sFSP_telegraph(parm,maxRna); %this new calcualtion had been added the 28.09.2021
    distRNA(distRNA<10e-20)=10e-20;%ATTENTION this is to remove extremely small probabilities
    if Bwidth>1
        %modified 05.01.2021
        %distRNAbinned=kron(eye((maxRna+1)/Bwidth),ones(1,Bwidth))*distRNA;%calculate the histogram of the distribution with bin width = Bwidth
        binnedmat=kron(eye(floor(length(distRNA)/Bwidth)),ones(1,Bwidth));
        if rem(length(distRNA),Bwidth)>0
            binnedmat(end+1,end+1:end+rem(length(distRNA),Bwidth))=ones(1,rem(length(distRNA),Bwidth));
        end
        distRNAbinned=binnedmat*distRNA;%calculate the histogram of the distribution with bin width = Bwidth
        LdistRNA=log(distRNAbinned); %log of the probability dist.
    else 
        LdistRNA=log(distRNA); %log of the probability dist.
    end
    LdistRNAtruncated=LdistRNA(1:length(dataHist)); %this is added to be able to sum up the pxlog(p) in the next line
    lddata=dataHist*LdistRNAtruncated;%sum the log probabilities weighted by the data histogram
    if isnan(lddata)==1 || isinf(lddata)==1
        like=10000000;
    else
        like=-sum(lddata);
    end
end