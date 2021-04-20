%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_fullSCR - 
%
%
%
% Other m-files required: fit_fct_all,fsp_ep.m, ll_ep.m, mean_rna.m,
% rate_matrix.m. fct_fitlongSCR.m
%
% Author: Gregory Roth
%
%   original version: 19.02.2020,
%   last version: 15.04.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------
%load the data
%-------------------------
load('fullSCRdata.mat')
load('bestfit_lowregime.mat')
parmL=x;

%-------------------------
%Fit
%-------------------------
%parameter for the fit
nbrun=10;%number of runs

nvec=2:1:11; %values for n+1

ltemp = [1.001;0;0;0;0;0;0;0;0;0];
utemp = [10;10;100;100;100;100;100;100;1000;3000]; 

%fit the full model
fitoption=1;
parmtofit=[0,1,1,1,0,0,1,1,0,1];%free parameters (n is fitted sequentially because it is an integer)
parmfixed=[nan,nan,nan,nan,parmL(1),parmL(2),nan,nan,parmL(3),nan];
for i=1:length(nvec)
    [x,fval,eflag,output,manymins]=fit_ep(fitoption,maxRnaL,dataValL,binsizeL,maxRnaH,dataValH,binsizeH,cpdata,rnadata,parmtofit,parmfixed,nvec(i),nbrun,ltemp,utemp);
    bestfit(i,:)=x;
end
