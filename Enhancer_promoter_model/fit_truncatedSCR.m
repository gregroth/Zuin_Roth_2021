%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fit_longSCR - 
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
load('truncatedSCRdata.mat')
load('bestfit_fullSCR.mat')


%-------------------------
%Fit
%-------------------------
%parameter for the fit
nbrun=100;%number of runs

ltemp = [1.001;0;0;0;0;0;0;0;0;0];
utemp = [10;10;100;100;100;100;100;100;1000;3000]; 

%fit the full model
fitoption=2;
parmtofit=[0,0,0,0,0,0,1,0,0,0];%choose the free parameters
parmfixed=bestfit_fullSCR;
parmfixed(parmtofit==1)=nan;
[x,fval,eflag,output,manymins]=fit_ep(fitoption,maxRnaL,dataValL,binsizeL,maxRnaH,dataValH,binsizeH,cpdata,rnadata,parmtofit,parmfixed,bestfit_fullSCR(1),nbrun,ltemp,utemp);
    
