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
function runPCA

global ica h_main st

status = 'Running PCA...';
set(h_main.settings.status.text, 'String', status), drawnow

% Finding dimensions of X:
ica.dimX = size(ica.X);

%% Mask data:

% Setup images in columns and make mask:
Xcolumns = reshape(ica.X, prod(ica.dimX(1:3)), ica.dimX(4));

% Find temporal variance components larger than mean  temporal variance:
maske = get(h_main.settings.pca.mask.method,'String');
maske = strtrim(maske(get(h_main.settings.pca.mask.method,'Value'),:));
switch maske
    case 'Variance'
        if size(ica.FileNames,1) > 20   % Limit for variance calculation due to memory reasons
            K = 20;
        else
            K = size(ica.FileNames,1);
        end
        varX = var(Xcolumns(:,1:K),0,2);

        ica.mask = find( varX > 2*mean(varX) );

    case 'Mean'
        ica.mask = find(mean(Xcolumns,2) > 3*mean(Xcolumns(:)));
        
    case 'User specified'
        ica.mask = spm_select(1,'mat','Choose file containing mask');
        
        if sum(size(ica.mask)) == 0
            error('Mask is not recognized')
        elseif ndims(ica.mask) == 3    % Mask is in 3D
            if numel(ica.mask) == size(Xcolumns,1)
                ica.mask = find(ica.mask); % Convert to indices
            else
                error(['Mask does not have the correct size (numel(mask) = '...
                    num2str(numel(ica.mask)) ' and numel(image) = '...
                    num2str(size(Xcolumns,1))])
            end
        elseif ndims(ica.mask) == 2     % Mask is correctly defined
        else
            error('Mask is not recognized')
        end
        
    case 'None'
        ica.mask = 1:size(Xcolumns,1);
end

if get(h_main.settings.pca.showmask,'Value')
    masken = zeros(size(Xcolumns,1),1);
    masken(ica.mask)=1;
    masken = reshape(masken,ica.dimX(1),ica.dimX(2),ica.dimX(3));
    temp = ica.Header(1);
    [PATHSTR,NAME,EXT] = FILEPARTS(temp.fname); 
    temp.fname = fullfile(PATHSTR, 'ica_img.img');
    masken = spm_write_vol(temp, masken);
    temp.st.fig = st.fig;
    st.fig = figure('units','normalized',...
        'Position',[0.26 0.05 .7 .9],...
        'name', 'Mask visualisation');
    spm_orthviews('Image',masken);
    colormap gray
    st.fig = temp.st.fig;
    st.h.orthviews = spm_orthviews('Image',ica.StructuralHeader,[.05 .45 .48 .48]);
end


ica.Xmasked = Xcolumns(ica.mask,:);

%% Highpass filter data
K.RT = str2double(ica.TR{1});               % Repetition time
K.row = 1:size(ica.Xmasked,2);  % Time points
K.HParam = 60;                  % Cut-off period [s]
K = spm_filter(K);              % Filter setup
X_T = ica.Xmasked';
X_T = spm_filter(K, X_T); % Filtering
ica.Xmasked = X_T';

%% Make PCA
% Subtract mean:
ica.Xmasked = ica.Xmasked - repmat(mean(ica.Xmasked,2), 1, ica.dimX(4));
% Perform SVD:
[ica.pca.U, ica.pca.S, ica.pca.V]= svd(ica.Xmasked, 0);

status = 'Done running PCA.';
set(h_main.settings.status.text, 'String', status)