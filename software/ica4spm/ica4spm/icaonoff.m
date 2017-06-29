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
function icaonoff

global h_main ica

if get(h_main.settings.pca.check,'Value')   % Perform PCA on
    set(h_main.settings.pca.mask.text,'Enable','on')
    set(h_main.settings.pca.mask.method,'Enable','on')
    set(h_main.settings.pca.showmask,'Enable','on')    
else                                        % Perform PCA off
    set(h_main.settings.pca.mask.text,'Enable','off')
    set(h_main.settings.pca.mask.method,'Enable','off')
    set(h_main.settings.pca.showmask,'Enable','off')        
end

if get(h_main.settings.ica.check,'Value')   % Perform ICA on
    set(h_main.settings.ica.algorithm.text,'Enable','on')
    set(h_main.settings.ica.algorithm.method,'Enable','on')
    set(h_main.settings.ica.tau.check,'Enable','on')
    set(h_main.settings.ica.tau.set,'Enable','off')
    set(h_main.settings.ica.type.text,'Enable','on')
    set(h_main.settings.ica.type.method,'Enable','on')
    set(h_main.settings.ica.bic.check,'Enable','on')
    set(h_main.settings.ica.ICs.text,'Enable','on')
    set(h_main.settings.ica.ICs.set,'Enable','on')
    if get(h_main.settings.ica.bic.check,'Value')   % Perform BIC on
        set(h_main.settings.ica.ICs.text,'String','Maximum number of ICs')
        set(h_main.settings.ica.ICs.set,'String','40')
        ica.ICs = str2double(get(h_main.settings.ica.ICs.set,'String'));
    else                                            % Perform BIC off
        set(h_main.settings.ica.ICs.text,'String','Specify number of ICs')
    end
else                                        % Perform ICA off
    set(h_main.settings.ica.algorithm.text,'Enable','off')
    set(h_main.settings.ica.algorithm.method,'Enable','off')
    set(h_main.settings.ica.tau.check,'Enable','off')
    set(h_main.settings.ica.tau.set,'Enable','off')
    set(h_main.settings.ica.type.text,'Enable','off')
    set(h_main.settings.ica.type.method,'Enable','off')
    set(h_main.settings.ica.bic.check,'Enable','off')
    set(h_main.settings.ica.ICs.text,'Enable','off')
    set(h_main.settings.ica.ICs.set,'Enable','off')
end