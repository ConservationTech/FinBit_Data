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

function icashow

global ica st h_main

status = 'Running ICA visualization...';
set(h_main.settings.status.text, 'String', status), drawnow

% Find invoking function:
a = dbstack;
if length(a) > 1
    comp = str2double(get(st.h.ic.component,'String'));
else
    if gcbo == st.h.ic.next && str2double(get(st.h.ic.component,'String')) < ica.ICs
        comp = str2double(get(st.h.ic.component,'String'))+1;
        set(st.h.ic.component,'String',num2str(comp))
    elseif (gcbo == st.h.ic.previous) && str2double(get(st.h.ic.component,'String')) > 1
        comp = str2double(get(st.h.ic.component,'String'))-1;
        set(st.h.ic.component,'String',num2str(comp))
    else
        comp = str2double(get(st.h.ic.component,'String'));
    end
end

% Determine fraction of blobs to show:
frac = [ica.ICsvisibility/100 1-ica.ICsvisibility/100];
voxels = prod(ica.dimX(1:3));
x = reshape(ica.image(:,:,:,comp),voxels,1);
[q,k] = sort(x);

ica.ica_min = zeros(ica.dimX(1:3));
ica.ica_max = ica.ica_min;

indx_low=k(1:floor(voxels*frac(1)));
indx_high=k(ceil(voxels*frac(2)): voxels);

ica.ica_min(indx_low)=1;
ica.ica_max(indx_high)=1;

temp = ica.StructuralHeader;
[PATHSTR,NAME,EXT] = FILEPARTS(temp.fname); 
temp.fname = fullfile(PATHSTR, 'ica_img_min.img');
ica.ica_min = spm_write_vol(temp, ica.ica_min);

temp.fname = fullfile(PATHSTR, 'ica_img_max.img');
ica.ica_max = spm_write_vol(temp, ica.ica_max);

spm_orthviews('RemoveBlobs',st.h.orthviews)
spm_orthviews('addcolouredimage',st.h.orthviews,ica.ica_min, [0 0 1])
spm_orthviews('addcolouredimage',st.h.orthviews,ica.ica_max, [1 0 0])
spm_orthviews('Redraw')

TR = str2double(ica.TR{1});

if isempty(ica.paradigm.paradigm)   % Plot current IC without paradigm
    if get(st.h.timeplot.showcurrent,'Value')
        plot(st.h.timeplot.axes, 1:TR:TR*ica.dimX(4), ica.S(comp,:),'b')
        set(st.h.timeplot.axes,'Xlim',[0 TR*ica.dimX(4)])
        legend(st.h.timeplot.axes, ['Comp. ' num2str(comp)])
        plot(st.h.freqplot.axes, linspace(0,1/TR/2,size(ica.componentfreq,2)), ica.componentfreq(comp,:))
        set(st.h.freqplot.axes,'Xlim',[0 1/TR/2])
    else                            % Plot all ICs without paradigm
        max_amp = max(ica.S(:));
        amp_shift = repmat(max_amp*(1:ica.ICs)',1,size(ica.S,2));
        plot(st.h.timeplot.axes, 1:TR:TR*ica.dimX(4), ica.S+amp_shift)
        set(st.h.timeplot.axes,'Xlim',[0 TR*ica.dimX(4)])
        legend(st.h.timeplot.axes, [repmat('Comp. ',ica.ICs,1) num2str((1:ica.ICs)')])
        plot(st.h.freqplot.axes, linspace(0,1/TR/2,size(ica.componentfreq,2)), ica.componentfreq)
        set(st.h.freqplot.axes,'Xlim',[0 1/TR/2])
    end
else                                % Plot current ICs with paradigm
    if get(st.h.timeplot.showcurrent,'Value')
        plot(st.h.timeplot.axes, 1:TR:TR*ica.dimX(4), ica.S(comp,:),'b',...
            ica.paradigm.t,ica.paradigm.paradigm,':r')
        set(st.h.timeplot.axes,'Xlim',[0 TR*ica.dimX(4)])
        legend(st.h.timeplot.axes, ['Comp. ' num2str(comp)], 'Paradigm')
        plot(st.h.freqplot.axes, linspace(0,1/TR/2,size(ica.componentfreq,2)), ica.componentfreq(comp,:))
        set(st.h.freqplot.axes,'Xlim',[0 1/TR/2])
    else                            % Plot all ICs with paradigm
        max_amp = max(ica.S(:));
        amp_shiftS = repmat(max_amp*(1:ica.ICs)',1,size(ica.S,2));
        amp_shiftp = repmat(max_amp*(1:ica.ICs)',1,size(ica.paradigm.paradigm,2));
        axes(st.h.timeplot.axes)
        plot(1:TR:TR*ica.dimX(4), ica.S + amp_shiftS), hold on
        plot(ica.paradigm.t, repmat(ica.paradigm.paradigm,ica.ICs,1) ...
            + amp_shiftp, ':r')
        set(st.h.timeplot.axes,'Xlim',[0 TR*ica.dimX(4)])
        legend(st.h.timeplot.axes, strvcat([repmat('Comp. ',ica.ICs,1) num2str((1:ica.ICs)')], ...
            'Paradigm'))
        hold off
        plot(st.h.freqplot.axes, linspace(0,1/TR/2,size(ica.componentfreq,2)), ica.componentfreq)
        set(st.h.freqplot.axes,'Xlim',[0 1/TR/2])
    end
end
title(st.h.timeplot.axes, 'Time series')
xlabel(st.h.timeplot.axes, 'Time [s]')
title(st.h.freqplot.axes, 'Component frequency content')
xlabel(st.h.freqplot.axes, 'Frequency [Hz]')

set(st.h.ic.variance.text, 'String',...
    sprintf(['Component variance is: ',...
    num2str(ica.variancecomponents(comp),2),...
    '\nTotal variance is: ' num2str(ica.totalvariance,2),...
    '\nRelative variance for component ' num2str(comp) ' is: ',...
    num2str(100*ica.variancecomponents(comp)/ica.totalvariance,2) ' %%']))

status = 'Done running ICA visualization...';
set(h_main.settings.status.text, 'String', status)
