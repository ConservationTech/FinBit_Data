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
function loadparadigm

global ica st

if get(st.h.timeplot.showparadigm, 'Value')
    i = strfind(ica.FileNames(1,:), '\');
    i = i(end);

    current = pwd;
    cd(ica.FileNames(1,1:i))
    files = dir('*.mat');
    cd(current)
    existparadigm = zeros(length(files),1);
    for n = 1:length(files)
        existparadigm(n) = strcmp(files(n).name,'paradigm.mat');
    end
    if sum(existparadigm)
        ica.paradigmpathname = [ ica.FileNames(1,1:i) files(find(existparadigm)).name];
        load(ica.paradigmpathname);
        ica.paradigm.paradigm = paradigm;
        ica.paradigm.t = t;
    else
        set(st.h.timeplot.showparadigm, 'Value',0)
        ica.paradigm.paradigm = [];
        ica.paradigm.t = [];
        errordlg(['There is no file called paradigm.mat in the directory ' ica.FileNames(1,1:i)])
        return
    end
else
    ica.paradigm.paradigm = [];
    ica.paradigm.t = [];
end

icashow