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

function tauonoff

global h_main

if get(h_main.settings.ica.algorithm.method,'Value') == 1 && ...
        get(h_main.settings.ica.tau.check,'Value') == 1
    set(h_main.settings.ica.tau.check,'enable','on')
    set(h_main.settings.ica.tau.set,'enable','on')
elseif get(h_main.settings.ica.algorithm.method,'Value') == 1 && ...
        get(h_main.settings.ica.tau.check,'Value') == 0
    set(h_main.settings.ica.tau.check,'enable','on')
    set(h_main.settings.ica.tau.set,'enable','off')
elseif get(h_main.settings.ica.algorithm.method,'Value') == 2
    set(h_main.settings.ica.tau.check,'tooltipstring','Only relevant with Molgedey-Schuster')
    set(h_main.settings.ica.tau.check,'enable','off')
    set(h_main.settings.ica.tau.set,'enable','off')
    set(h_main.settings.ica.tau.set,'string','')
end
