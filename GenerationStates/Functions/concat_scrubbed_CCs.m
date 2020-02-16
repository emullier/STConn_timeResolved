

function varargout = BD_concat_scrubbed_CCs(path_fmriprep, path_stconn, scale, varargin)

% 
% Syntax :
%  varargout = concat_scrubbed_CCs(path_fmriprep, path_stconn, scale, varagin);
%
% This script concatenates the CCs of a group of subjects from the output results of the spatio-temporal connectome [Griffa2017].
% It also scrubs them to remove the small CCs and the CCs with time points
% affected by motion (according to FD and DVARS scrubbing parameters).
%
% Input Parameters:
%      scale           :  scale number of Lausanne2008 parcellation (1,2,3,4 or 5)
%      path_stconn     :  path to output of the spatio-temporal connectome 
%      path_fmriprep   :  path to output of the fmri preprocessing 
%      thFD            :  threshold for FD for the scrubbing of CCs affected by motion
%      thDVARS         :  threshold for DVARS for the scrubbing of CCs affected by motion
%      min_height      :  minimal height of the CCs
%      min_width       :  minimal width of the CCs
%      
%
% Output Parameters:
%     list_CCs         :  concatenated list of CCs for all the subjects
%     list_subs        :  list of subjects number corresponding to the list of CCs
%
% Related references:
%
%  Usage: 
%   [list_CCs, list_subs] = concat_scrubbed_CCs('/media/localadmin/Windows/Emeline/STConn_SCHZ/fmri_preprocessing', '/media/localadmin/Windows/Emeline/STConn_SCHZ/STConn_scale4', 4)
%   [list_CCs, list_subs] = concat_scrubbed_CCs('/media/localadmin/Windows/Emeline/STConn_SCHZ/fmri_preprocessing', '/media/localadmin/Windows/Emeline/STConn_SCHZ/STConn_scale4', 4, 'min_height', 6)
%  
% See also:
%   CCs2BrainMaps.m 
%
%__________________________________________________
% Authors: Emeline Mullier
% Connectomics Lab, Department of Radiology, Lausanne University Hospital
% April 4th 2019
% Version $3.0



%% ====================== Checking input parameters ===================== %
if nargin<3 % the indispensable input arguments are not provided
    error('Three inputs are mandatory: path_fmriprep, path_stconn and scale');
else
    thFD       = 0.5;
    thDVARS    = 20;
    min_height = 6;
    min_width  = 2;
end

% deal with the input arguments
if nargin<3 % the indispensable input arguments are not provided
    error('Three inputs are mandatory: path_fmriprep, path_stconn and scale');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'thFD' % FD threshold for CCs scrubbing
                    thFD=varargin{2};
                case 'thDVARS' % DVARS threshold for CCs scrubbing
                    thDVARS=varargin{2};
                case 'min_height' % minimal CC height
                    min_height=varargin{2};
                case 'min_width' % minimal CC width
                    min_width=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end
%% ================= End of Checking input parameters =================== %


list_CCs = [];
list_subs = [];
list_nb_CCs = [];
list_FD = [];
list_DVARS = [];

list_sub = dir(fullfile(path_stconn, '*sub*'));

clear subs_ids;
for i=1:length(list_sub)
    
    path_ses = sprintf('%s/%s', path_stconn, list_sub(i).name);
    path_CC = sprintf('%s/CC.mat', path_ses);
    path_ts = sprintf('%s/ts_zscore.mat', path_ses);
   
    path_FD    = sprintf('%s/%s/FD.mat', path_fmriprep, list_sub(i).name);
    path_DVARS = sprintf('%s/%s/DVARS.mat', path_fmriprep, list_sub(i).name);

    if exist(path_CC)
        subs_ids{i} = list_sub(i).name;
        load(path_CC); load(path_ts); 
        load(path_FD); load(path_DVARS);
        for j=1:length(CC)
            tp_ind = CC{j}.tp; tp_ind(find(tp_ind==size(ts,2)))=[];
            if (CC{j}.height > min_height) && (CC{j}.width > min_width)
                if (max(FD(unique(tp_ind)))<thFD) && (max(DVARS(unique(tp_ind)))<thDVARS)
                    list_CCs  = [list_CCs CC{j}];
                    list_subs = [list_subs; i];
                    list_nb_CCs = [list_nb_CCs; j];
                    list_FD = [list_FD; FD];
                    list_DVARS = [list_DVARS; DVARS];

                end
            end
        end       
    else
        error(sprintf('The CCs do not exist for %s', list_sub(i).name));
    end
end



%%  ------------------ Output parameters -------------------------------- %
varargout{1} = list_CCs;     % list of CCs
varargout{2} = list_subs;    % list of corresponding subjects
varargout{3} = list_nb_CCs;  % CC corresponding to the time point number
varargout{4} = list_FD;      % CC corresponding to the time point number
varargout{5} = list_DVARS;   % CC corresponding to the time point number



end