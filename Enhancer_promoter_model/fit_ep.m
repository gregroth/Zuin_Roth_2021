

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_ep - This function fits the enhancer-promoter model with specific
%number of regulatory steps to a specific set of data
%
%
%Input
%
%           - countL
%           - countH
%           - cpdata
%           - rnadata,
%           - cvdata,
%           - cpalldata,
%           - cpFISH
%           - cvFISH,
%           - parmvec
%           - nbrun,
%           - parmfixed
%           - nbs
%
%Output
%
%           -  dist.join= vector containing the distribution of the state
%              of the promoter, the state of the enhancer, and the number of RNA
%           -  dist.rna= vector containin
%
% Other m-files required: fit_fct_all,fsp_ep.m, ll_ep.m, mean_rna.m,
% rate_matrix.m. fct_fitlongSCR.m
%
% Author: Gregory Roth
%
%   original version: 19.02.2020,
%   last version: 15.04.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,fval,eflag,output,manymins]=fit_ep(fitoption,maxRnaL,dataValL,binsizeL,maxRnaH,dataValH,binsizeH,cpdata,rnadata,parmvec,parmfixed,nbs,nbrun,ltemp,utemp)

N=length(rnadata); %number of data points
sigma=1; %standard deviation for the mean rna model

x0temp = 9*rand(1,length(parmvec))+.00001; %initial parameters
x0temp(1)=x0temp(1)+1;
x0temp(2)=.01;
    
%-------------------------
%Fit
%-------------------------
t=10;

%define the model for the mean
if parmvec==[0,1,1,1,0,0,1,1,0,1]
    parmL=parmfixed;
    parmL(isnan(parmL))=[];
    rna_model=@(P,cp) mean_rna([nbs,P(1:3),parmL(1:2),P(4:5),parmL(3),P(6)],cp);
elseif parmvec==[0,1,0,0,0,0,0,0,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),P,parmfixed(3),parmfixed(4),parmfixed(5),parmfixed(6),parmfixed(7),parmfixed(8),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,1,1,0,0,0,0,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),P(1),P(2),parmfixed(5),parmfixed(6),parmfixed(7),parmfixed(8),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,1,0,0,0,0,0,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),P(1),parmfixed(4),parmfixed(5),parmfixed(6),parmfixed(7),parmfixed(8),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,0,1,0,0,0,0,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),P(1),parmfixed(5),parmfixed(6),parmfixed(7),parmfixed(8),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,0,0,0,0,1,1,0,1]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),parmfixed(4),parmfixed(5),parmfixed(6),P(1),P(2),parmfixed(9),P(3)],cp);
elseif parmvec==[0,0,1,1,0,0,1,1,0,1]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),P(1),P(2),parmfixed(5),parmfixed(6),P(3),P(4),parmfixed(9),P(5)],cp);
elseif parmvec==[0,0,1,1,0,0,1,1,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),P(1),P(2),parmfixed(5),parmfixed(6),P(3),P(4),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,1,0,0,0,1,1,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),P(1),parmfixed(4),parmfixed(5),parmfixed(6),P(2),P(3),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,0,1,0,0,1,1,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),P(1),parmfixed(5),parmfixed(6),P(2),P(3),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,1,0,0,0,1,1,0,1]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),P(1),parmfixed(4),parmfixed(5),parmfixed(6),P(2),P(3),parmfixed(9),P(4)],cp);
elseif parmvec==[0,0,0,1,0,0,1,1,0,1]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),P(1),parmfixed(5),parmfixed(6),P(2),P(3),parmfixed(9),P(4)],cp);
elseif parmvec==[0,0,0,0,0,0,1,1,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),parmfixed(4),parmfixed(5),parmfixed(6),P(1),P(2),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,0,0,0,0,1,0,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),parmfixed(4),parmfixed(5),parmfixed(6),P(1),parmfixed(8),parmfixed(9),parmfixed(10)],cp);
elseif parmvec==[0,0,0,0,0,0,0,0,0,1]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),parmfixed(4),parmfixed(5),parmfixed(6),parmfixed(7),parmfixed(8),parmfixed(9),P(1)],cp);
elseif parmvec==[0,0,0,0,0,0,0,1,0,0]
    rna_model=@(P,cp) mean_rna([ceil(parmfixed(1)),parmfixed(2),parmfixed(3),parmfixed(4),parmfixed(5),parmfixed(6),parmfixed(7),P(1),parmfixed(9),parmfixed(10)],cp);
end    




x0 = x0temp(parmvec==1);
l = ltemp(parmvec==1);
u = utemp(parmvec==1);

%calculate the total log likelihood function
if fitoption==1
    LL_mean = @(x)N*log(sigma*sqrt(2*pi))+sum((rnadata-rna_model(x,cpdata)).^2)/(2*sigma^2);
    LL_high = @(x)ll_ep([nbs,x(1:3),parmL(1:2),x(4:5),parmL(3),x(6)],1,maxRnaH,dataValH,binsizeH,t); 
    LL_tot=@(x)LL_mean(x)+LL_high(x);
elseif fitoption==2
    LL_mean = @(x)N*log(sigma*sqrt(2*pi))+sum((rnadata-rna_model(x,cpdata)).^2)/(2*sigma^2);
    LL_tot=@(x)LL_mean(x);
elseif fitoption==3
    LL_low = @(x)ll_low(x,maxRnaL,t,dataValL,binsizeL); 
    LL_tot=@(x)LL_low(x);
end

%set the global optimzation problem
    problem = createOptimProblem('fmincon',...
     'objective',@(x)LL_tot(x),...
     'x0',x0,'lb',l,'ub',u);

%run the optimization
ms = MultiStart('PlotFcns',@gsplotbestf);
    
[x,fval,eflag,output,manymins] = run(ms,problem,nbrun);
end
