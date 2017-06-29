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


function xyz_pos

global st

if ~isempty(gcbo)
    if gcbo == st.h.pos.horizontal || gcbo == st.h.pos.coronal || gcbo == st.h.pos.sagittal
        centre = [str2double(get(st.h.pos.sagittal,'String')),...
            str2double(get(st.h.pos.coronal,'String')),...
            str2double(get(st.h.pos.horizontal,'String'))];
        spm_orthviews('Reposition',centre)
    else
        set(st.h.pos.coronal,'String',num2str(st.centre(2),3))
        set(st.h.pos.sagittal,'String',num2str(st.centre(1),3))
        set(st.h.pos.horizontal,'String',num2str(st.centre(3),3))
    end
end
