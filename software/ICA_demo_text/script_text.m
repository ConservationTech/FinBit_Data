% script to classify text using ICA
% by Thomas Kolenda DTU,IMM 2000,2002 version 2

close all
clear all
format compact

% settings
ClassFrac=0;               % Reject frac for classification [0..1[, zero means no reject class
keywordstreshold=.30;      % Treshold for finding keywords ]0..1[

% load data
load med_5class            % File holds term/document matrix X, terms, target labels
X=termDocNorm;
Targets=cell2mat(base.qrels);

% PCA
[T,L,D]=svd(X,0);
DL=L*D';

% Plot first five PC components against each other with given labels
figure(1); plotComp(DL,base,'PC comps with labels');    

% Estimating number of components by BIC
Xb.U=T; Xb.S=L; Xb.V=D;
BIC_P=icaML_bic(X,1:10);
[maxBIC_P,Components]=max(BIC_P);

disp(sprintf('Estimated %d components using BIC with probability %0.2f\n',Components,maxBIC_P));

% ICA decomp
[S,A]=icaML(DL(1:Components,:));

% Plot first five IC components against each other with given labels as color
figure(2); plotComp(S,base,'IC comps with labels');    

% Find positive direction, Find soft components and classify by angle
[S,A]=flipcomp_ica(S,A);

% Use if outputs should be regarted as probabilities P(S|K)
%[S,A]=softmax_ica(S,A);    

% Find classes
Estimats=classifyer_angel_ica(S,ClassFrac);

% Plot first five soft and fliped IC components against each other with estimated class colors
figure(3); plotComp(S,Estimats,'IC comps with estimated classes');    

% Plot BIC estimat
figure(4); bar(1:10,BIC_P); ylabel('P(K)'); xlabel('K'); title('BIC estimat')

% Cals keywords and build confusion matrix
keystr=calckeywords(T,A,terms,keywordstreshold);
Confusion_Matrix=confusionMatrix(Targets,Estimats)

