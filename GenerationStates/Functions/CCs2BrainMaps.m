function varargout = BD_CCs2BrainMaps(list_CCs, list_subs, scale, nb_tp, nROIs)

% 
% Syntax :
%  varargout = BD_CCs2BrainMaps(list_CCs, list_subs, scale);
%
% This script converts a list of CCs into a series of brain
% configurations, activations (activated ROIs) maps 
%
% Input Parameters:
%      list_CCs         :  concatenated list of CCs for all the subjects
%      list_subs        :  list of subjects number corresponding to the list of CCs
%      scale            :  scale number of Lausanne2008 parcellation (1,2,3,4 or 5)
%      nb_tp            :  number of frames of the acquisition.
%      nROIs            :  number of regions in the parcellation
%      
% Output Parameters:
%      brain_maps          :  the series of brain maps
%
%
% Related references:
%
%  Usage: 
%   [brain_maps] = BD_CCs2BrainMaps(list_CCs, list_subs, scale);
%  
% See also:
%   BD_concat_scrubbed_CCs.m, BD_StatesIdentification.m 
%
%__________________________________________________
% Authors: Emeline Mullier
% Connectomics Lab, Department of Radiology, Lausanne University Hospital
% September 9th 2019
% Version $3.0



%% ====================== Checking input parameters ===================== %
if nargin<3 % the indispensable input arguments are not provided
    error('Three inputs are mandatory: list_CCs, list_subs and scale');
end
%% ================= End of Checking input parameters =================== %


nCCs = length(list_CCs);
ntp_concat = nb_tp * max(list_subs);
  
brain_maps = zeros(nROIs, ntp_concat);

for i=1:nCCs
    idxs = find(list_CCs(i).tp>nb_tp);
    ls_nodes = list_CCs(i).nodes(idxs);
    for j=1:size(list_CCs(i).edges,1)
        if ~ismember(list_CCs(i).edges(j,1),ls_nodes) && ~ismember(list_CCs(i).edges(j,2),ls_nodes)
            [x1,y1] = ind2sub([nROIs, nb_tp], list_CCs(i).edges(j,1));
            [x2,y2] = ind2sub([nROIs, nb_tp], list_CCs(i).edges(j,2));
            if (y1<=nb_tp) && (y2<=nb_tp)
            %if y1==y2
                brain_maps(x1,y1+nb_tp*(list_subs(i)-1)) = 1;
                brain_maps(x2,y2+nb_tp*(list_subs(i)-1)) = 1;
            %end
            end
        end
    end
end


%%  ------------------ Output parameters -------------------------------- %
varargout{1} = brain_maps;         % series of brain maps corresponding to a list of concatenated CCs (across subjects)


end