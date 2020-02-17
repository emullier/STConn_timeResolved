


%%% =========================================================================================================================
%%%
%%% This script allows to compare the states obtained with
%%% 'GenerationSates' to the Yeo's 7 functional systems [Yeo2011].
%%% 
%%% Emeline Mullier - February 2020
%%% Connectomics Lab, Department of Radiology, Lausanne University Hospital
%%%
%%% Manuscript
%%% 'Functional brain dynamic are shaped by connectome n-simplicial organization'
%%% Emeline Mullier, Jakub Vohryzek, Alessandra Griffa, Yasser Alemàn-Gómez, Celia Hacker, Kathryn Hess, Patric Hagmann
%%%
%%% ==========================================================================================================================


addpath(genpath('BCT'));
addpath(genpath('Plot_Brains_Kit'));

cd '/home/localadmin/Documents/STConn_Analyses/STConn_timeResolved';
working_dir = 'CompStates';
addpath(genpath(working_dir));


%%%================== %%%
%%% Parameters        %%%
%%%================== %%%

load('labels_index_CORTICAL_Laus2008_all_scales.mat'); % indices of the cortical regions
scale = 4;   % scale of the parcellation
nb_tp = 275; % number of time points/ TR
nROIs = 463; % number of regions
thNoise = 2; % threshold to consider a community as noise
fslabel = 12; fstitle = 16; % fontsize for plotting

export = 1; % choose if the figures need to be exported or not

fig_path = 'Figures/'; % path to save the figures
suffixe = 'STConn_g1'; % suffixe used to save the results
load('GenerationStates/Data/CCsPrep_CCs.mat'); load('Data/Louvain_g1_init20.mat'); % load the analysis

%%% Removing noisy clusters
idxsNoise = find(accumarray(CIU, ones(size(CIU)))>thNoise);
mat_clusters = mat_clusters(ixc{scale},idxsNoise);


%%%================%%%
%%% YEO SIMILARITY %%%
%%%================%%%

%%% Load Yeo-Lausanne2008 correspondance
labels_yeo = {'VIS', 'SM', 'DA', 'VA', 'LIMB', 'FP', 'DMN'};
nbYeo = length(labels_yeo);
load(sprintf('%s/Functions/Lausanne2Yeo.mat', working_dir));
Lausanne2Yeo = Lausanne2Yeoh.scale4;
Lausanne2Yeo = Lausanne2Yeo(:,2:end);

%%% Correlation betwen Yeo's systems and the states
MatSimYeo = zeros(size(mat_clusters,2), nbYeo);
for i=1:size(mat_clusters,2)
    for j=1:nbYeo
        YeoSys = Lausanne2Yeo(:,j);
        r = corrcoef(mat_clusters(:,i), YeoSys);
        MatSimYeo(i,j) = r(1,2);
    end
end
figure; imagesc(MatSimYeo);  axis square; colorbar('location', 'southoutside');
ylabel(sprintf('%s',suffixe), 'fontsize', fslabel, 'interpreter', 'latex');
set(gca, 'XTick', [1:nbYeo], 'XTickLabel', labels_yeo, 'fontsize', fslabel, 'TickLabelInterpreter', 'latex');
xtickangle(45); title('Correlation', 'fontsize', fstitle, 'interpreter', 'latex'); colorbar;
Figurename = [working_dir '/' fig_path sprintf('MatSim_%s_Yeo.csv', suffixe)];
csvwrite(Figurename, MatSimYeo');
if export
    Figurename = [working_dir '/' fig_path sprintf('MatSim_%s_Yeo.tiff', suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300');
    close;
end




