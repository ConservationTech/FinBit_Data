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

function runICA

global h_main ica

status = 'Running ICA...';
set(h_main.settings.status.text, 'String', status), drawnow

% Represent the data matrix in terms of K principal components
Xreduced = ica.pca.U(:, 1:ica.ICs)' * ica.Xmasked;

draw = 0;

switch ica.Algorithm
    case 'Molgedey-Schuster '
        if get(h_main.settings.ica.tau.check,'Value') == 1  %
            if isempty(str2double(get(h_main.settings.ica.tau.set,'string')))
                ica.tau = 0;
            else
                ica.tau = str2double(get(h_main.settings.ica.tau.set,'string'));
            end
        else
            ica.tau = 0;
        end

        switch ica.Type
            case 'Spatial '
                [ica.S,ica.A,ica.ll,ica.tau]=icaMS(Xreduced',ica.tau,draw);
                ica.S = ica.S';
                ica.A = ica.A';
            case 'Temporal'
                [ica.S,ica.A,ica.ll,ica.tau]=icaMS(Xreduced,ica.tau,draw);
        end

    case 'Maximum Likelihood'
        switch ica.Type
            case 'Spatial '
                [ica.S,ica.A,ica.U,ica.ll,ica.MLinfo]=icaML(Xreduced,0);
                ica.S = ica.S';
                ica.A = ica.A';
            case 'Temporal'
                [ica.S,ica.A,ica.U,ica.ll,ica.MLinfo]=icaML(Xreduced,0);

            otherwise
                error('The method is not implemented yet')
        end
        
    case 'FastICA           '
        [ica.S ica.A W] = fastica(Xreduced);
    
    case 'Jade              '
        ica.B = MatlabjadeR(Xreduced);
        ica.A = ica.B;  % Not correct. What to do???
        ica.S = ica.B*Xreduced;
end

ica.image = zeros(ica.dimX(1),ica.dimX(2),ica.dimX(3),ica.ICs);
ica.variancecomponents = zeros(ica.ICs,1);
for comp = 1:ica.ICs
    %     ica.Energy(j)=sum(ica.A(:,comp).*ica.A(:,comp));
    x = ica.pca.U(:,1:ica.ICs)*ica.A(1:ica.ICs, comp);
    vol = zeros(prod(ica.dimX(1:3)), 1);
    vol(ica.mask) = x;
    ica.image(:,:,:,comp) = reshape(vol, ica.dimX(1), ica.dimX(2), ica.dimX(3));
    ica.variancecomponents(comp) = 1/(prod(ica.dimX)) * sum(ica.A(:, comp).^2) * sum(ica.S(comp, :).^2);
    
    ica.componentfreq(comp,:) = pwelch(ica.S(comp,:));
end

ica.totalvariance = var(ica.image(ica.mask));

status = 'Done running ICA...';
set(h_main.settings.status.text, 'String', status)
