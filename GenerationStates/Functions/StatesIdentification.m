




function varargout = StatesIdentification(brain_maps, varargin)

% 
% Syntax :
%  varargout = StatesIdentification(brain_maps, varargin)
%
% This script computes a Louvain community detection algorithm on pairwise
% similarity matrix from a series of brain maps. The similarity matrix is
% based on the cosine distance between the brain maps. The function
% computes also the centroids of each community and the attribution of each
% brain map to its respective centroid. 
%
% Input Parameters:
%      brain_maps      :  a series of brain maps
%      nTRs            :  number of time points per subject
%      nb_init         :  number of initialization for the Louvain community detection algorithm
%      gamma           :  resolution parameter for the community detection
%      
%
% Output Parameters:
%     CIU              : Community attribution of each non null brain map
%     mat_clusters     : Centroids (average brain map)of the communities
%     mat_Q            : Modularity for each iteration
%     
%     
%
% Related references:
%
%  Usage: 
%   [CIU, mat_clusters, mat_Q] = StatesIdentification(brain_maps)
%   [CIU, mat_clusters, mat_Q] = StatesIdentification(brain_maps, 'nTRs', 275. 'nb_init', 5, 'gamma', 1.5)
%  
% See also:
%   concat_scrubbed_CCs.m, BD_CCsNoOverlaps2BrainMaps.m
%
%__________________________________________________
% Authors: Emeline Mullier
% Connectomics Lab, Department of Radiology, Lausanne University Hospital
% April 4th 2019
% Version $3.0


%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory: brain_maps');
else
    nTRs = 275;
    nb_init = 20;
    gamma = 1;
end

% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory: brain_maps');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'nTRs' % number of time points of the fMRI series for each subject
                    nTRs=varargin{2};
                case 'nb_init' % number of initialization for Louvain community detection
                    nb_init=varargin{2};
                case 'gamma' % resolution parameter for Louvain community detection
                    gamma=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% ================= End of Checking input parameters =================== %


%%% Compute the similarity matrix between all the non null brain maps
brain_maps_nonzeros = brain_maps(:,sum(brain_maps,1) ~= 0);
sim_cosine_matrix = pdist([brain_maps_nonzeros'],'cosine'); 
sim_mat = squareform(sim_cosine_matrix); sim_mat = ones(size(sim_mat)) - sim_mat;
sim_mat = sim_mat - eye(size(sim_mat));
fprintf('Similarity matrix computed \n');

%%% Louvain Community Detection
partitionMatrix = zeros(length(sim_mat), nb_init);
mat_Q = zeros(1, nb_init);
D = zeros(length(sim_mat), length(sim_mat));
CIU = zeros(length(sim_mat));
for i=1:nb_init
    [indComm,Q] = community_louvain(sim_mat,gamma); % get community assignments
    mat_Q(1,i) = Q;
    partitionMatrix(:,i) = indComm;
end
D = agreement(partitionMatrix);
CIU = consensus_und(D,0,5);
fprintf('Community attribution done \n');


%%% Average brain maps: States
mat_clusters = [];
for k=1:max(CIU)
    cluster = mean(brain_maps_nonzeros(:,find(CIU==k)),2);
    mat_clusters = [mat_clusters cluster];
end
fprintf('Average brain maps generated \n');


%%  ------------------ Output parameters -------------------------------- %
varargout{1} = CIU;          % Community attribution of the brain map series
varargout{2} = mat_clusters; % Centroids of the communities
varargout{3} = mat_Q;        % Modularity of the different iterations


end