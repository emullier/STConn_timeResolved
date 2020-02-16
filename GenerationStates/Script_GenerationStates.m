


%%% =========================================================================================================================
%%%
%%% This script allows to generate the states from the
%%% structurally-constrained time resolved brain maps generated from the
%%% spatio-temporal connectome [Griffa2017][Mullier,Submitted] on the
%%% provided 'DataExample' composed of 3 healthy subjects processed with
%%% the STConn [Griffa2017]
%%% 
%%% Emeline Mullier - December 2019
%%% Connectomics Lab, Department of Radiology, Lausanne University Hospital
%%%
%%% Manuscript
%%% 'The role of structural connectivity in functional brain dynamics' 
%%% Emeline Mullier, Jakub Vohryzek, Alessandra Griffa, Yasser Alemàn-Gómez, Celia Hacker, Kathryn Hess, Patric Hagmann
%%%
%%% ==========================================================================================================================


%addpath(genpath('/home/localadmin/Documents/STConn_Analyses/Plot_Brains_Kit'));


%% Set paths and paramaters
addpath(genpath('BCT')); %%% add path to the brain connectivity toolbox [Rubinov2010]
addpath(genpath(sprintf('%s',pwd)));

working_dir = 'GenerationStates';  % working directory
scale = 4;                         % scale of the parcellation Lausanne2008 [Cammoun2012]
nb_tp = 275;                       % number of fMRI volumes per subject
nROIs = 463;                       % number of regions in the parcellation

nb_init = 20;                      % number of iterations of the Louvain community detection
gamma = 1;                         % resolution of the Louvain community detection

CCsprep = 1;                       % boolean, to concatenate and scrub the CCs (according to fMRI motion parameters, FD and DVARS) and generate the matrix of structurally constrained time-resolved brain maps
Louvain = 1;                       % boolean, run the Louvain community detection
PlotStates = 1;                    % boolean, plot the average brain maps (=states)



%% ========== CCs PREPARATION ========== %%
preproc_path = 'DataExample/fMRI_preprocessing';
stconn_path = 'DataExample/STConn';
if CCsprep
    [list_CCs, list_subs] = concat_scrubbed_CCs(preproc_path, stconn_path, scale);
    brain_maps = CCs2BrainMaps(list_CCs, list_subs, scale, nb_tp, nROIs);
    save(sprintf('%s/Data/CCsPrep_CCs.mat',working_dir), 'list_CCs', 'list_subs', 'brain_maps');
else
    load(sprintf('%s/Data/CCsPrep_CCs.mat',working_dir))
end


%% ========== LOUVAIN COMMUNITY DETECTION ========== %%
if Louvain
        fprintf(sprintf('Gamma = %d \n', gamma));
        [CIU, mat_clusters, mat_Q] = StatesIdentification(brain_maps, 'nb_init', nb_init, 'gamma', gamma);
        save(sprintf('%s/Data/Louvain_g%d_init%d.mat',working_dir,gamma,nb_init), 'nb_init', 'gamma', 'CIU', 'mat_clusters', 'mat_Q');
end

%% ========== PLOT STATES ========== %%
if PlotStates
    if Louvain == 0
        load(sprintf('%s/Data/Louvain_g%d_init%d.mat', working_dir,gamma,nb_init));
    end
    load('PlotKit/labels_index_CORTICAL_Laus2008_all_scales.mat'); % indexes of the cortical regions for Lausanne2008 5 scales
    Surf_lh = Read_Surface('PlotKit/fsaverage/surf/lh.pial');
    Surf_rh = Read_Surface('PlotKit/fsaverage/surf/rh.pial');
    annotr = sprintf('PlotKit/fsaverage/label/rh.lausanne2008.scale%d.annot', scale);
    annotl = sprintf('PlotKit/fsaverage/label/lh.lausanne2008.scale%d.annot', scale);
    Noise = accumarray(CIU, ones(size(CIU)));
    mat_clusters = mat_clusters(ixc{scale},find(Noise>50));
    Plot_Brains_1view(Surf_rh, Surf_lh, annotr, annotl, mat_clusters','colMap','jet');
end



