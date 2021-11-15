%% routine to fit apparent two-state model to facs and fish data 
%
%
% Other m-files required: fit_fct.m/ mean_rna_mixedparam.m/
% MLdist_2state.m/ MLtot_apparent.m/ sFSP_telegraph.m
%
% Author: Gregory Roth
%
%   original version: 19.02.2021,
%   last version: 19.09.2021%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for the fit
%number of runs
nbrun=10;
%free parameters (N,konB,konE,koff,mu,beta)
freeparam=[1,1,1,1,1,1];
fixedparam=[];
%selected clones for the fit, belongs to 1,2,...,7
selectedclones=[1,1,1,1,1,0,1];
%initial point
initialpt=[0.0427    3.6310   10.0483  400  0.0020    2.7784];
%% multistart fit
[x,fval,eflag,output,manymins]=fit_fct(nbrun,freeparam,fixedparam,selectedclones,initialpt);
%extract the bestfit parameters orderd by goodness
Xtemp=arrayfun(@(x)x.X,manymins,'UniformOutput',false);
Xlist=cell2mat(Xtemp);
Xmat=reshape(Xlist,sum(freeparam),length(manymins))';
bestfit_param=zeros(1,6);
bestfit_param(freeparam==1)=x;
bestfit_param(freeparam==0)=fixedparam;

