%% fit_fct -- This function fits the variable apparent two-state model to facs and fish data
%
% Inputs:
%    nbrun - number of runs
%    freeparam - vector of size 6 with 1 for free parameter and 0 elsewhere
%    fixedparam - vector with the value of the fixed parameters
%    selectedclones - vector of size 7 with 1 if the clone is selected for
%                     the fit and 0 if it is not
%    initialpt - starting point for the optimization
%
% Other m-files required:  mean_rna_mixedparam.m/ ML_telegraph_hist.m/ MLdist.m/ sFSP_telegraph.m
% MAT-files required: fish_data.mat/ dataToFit.mat
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.02.2021%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,fval,eflag,output,manymins]=fit_fct(nbrun,freeparam,fixedparam,selectedclones,initialpt)
%% data to fit
load('../data/dataToFit.mat')
%remove the non selected clones
dataVal(selectedclones==0)=[];
binsize(selectedclones==0)=[];
maxRna(selectedclones==0)=[];
meanfish(selectedclones==0)=[];
cpfish(selectedclones==0)=[];

%attention divide by 2 the maxRna
maxRna=ceil(maxRna./1.6);

%constrain bounds 
lbtemp = [0;0;0;0;0;0];
ubtemp = [10;100;1000;10000;100000;100];
lb = lbtemp(freeparam==1);
ub = ubtemp(freeparam==1);

%likelihood functions
MLF = @(x)MLtot_app(freeparam,fixedparam,x,maxRna,dataVal,binsize,trdata_l,cpdata_l,cpfish);%fish dist


problem = createOptimProblem('fmincon',...
     'objective',@(x)MLF(x),...
     'x0',initialpt,'lb',lb,'ub',ub);
 
ms = MultiStart('UseParallel',false,'Display','iter');
%ms = MultiStart('PlotFcns',@gsplotbestf);
%ms = MultiStart('UseParallel',false,'Display','iter');
%set up parallel processing
%parpool
 % euler = parcluster('local');
  %pool = parpool(euler,24);
[x,fval,eflag,output,manymins] = run(ms,problem,nbrun);

end
