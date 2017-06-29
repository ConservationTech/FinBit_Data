% =========================================================================
% ICA4SPM toolbox
%
% - IMM, Technical University of Denmark
%
% - version 1.0 
% T. Bjerre, J. Henriksen†, C.H. Nielsen, P.M. Rasmussen
%                    L.K. Hansen, K.H. Madsen
%
% Bibtex reference:
%
%  @inproceedings{biosignals2009,
%    author = {Bjerre, T. and Henriksen, J. and Rasmussen, P.M. and Nielsen, C.H. and Hansen, L.K. and Madsen, K.H },
%    title = {Unified ICA-SPM analysis of f{MRI} experiments - Implementation of an {ICA} graphical user interface for the {SPM} pipeline},
%    booktitle = {Proceedings of BIOSTEC - BIOSIGNALS 2009 Conference},
%    year = {2009},}
%
%
%
% =========================================================================

function runBIC

global ica h_main st

status = 'Running BIC...';
set(h_main.settings.status.text, 'String', status), drawnow

M = size(ica.pca.U,1);
N = size(ica.pca.V,1);
K = 1:ica.ICs;

switch ica.Algorithm
    case 'Molgedey-Schuster '
        [ica.bic.P,ica.bic.logP]=icaMS_bic(ica.pca,K,0)
    case 'Maximum Likelihood'
        [ica.bic.P,ica.bic.logP]=icaML_bic(ica.pca,K,0)        
end

[maximum, bicsuggestion] = max(ica.bic.logP);

figure
bar(K,ica.bic.logP,1,'FaceColor',[.8 .8 .8]); 
hold on
bar(bicsuggestion,maximum,1,'r')
akse = axis;
axis([akse(1) K(end)+.5 akse(3) akse(4)])
title('BIC model selection');
ylabel('log of P(\it{X}|\it{K})');
xlabel('Number of independent components \it{K}');

acceptbic = questdlg(sprintf(['Suggested number of ICs: %d\n',...
    'Do you wish to use this number of ICs?'], bicsuggestion),...
    'Accept BIC suggestion','Yes','No, specify other','default');

switch acceptbic
    case 'Yes'
        ica.ICs = bicsuggestion;
    otherwise
        ica.ICs = str2double(inputdlg('Specify number of ICs'));
        if ica.ICs > N
            ica.ICs = N;
            warndlg(sprintf(['Specified number exceeded possible number of ICs.\n'...
                'ICs is set to %d'], ica.ICs))
        end
end
set(st.h.ic.numberICs,'String',['Number of ICs: ' num2str(ica.ICs)])

status = 'Done running BIC...';
set(h_main.settings.status.text, 'String', status)
