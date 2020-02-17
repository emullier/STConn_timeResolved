


%%% =========================================================================================================================
%%%
%%% This script allows to dynamic behaviour measures of the states with the
%%% provided 'DataExample' composed of 3 healthy subjects processed with
%%% the STConn [Griffa2017]
%%% 
%%% Emeline Mullier - December 2019
%%% Connectomics Lab, Department of Radiology, Lausanne University Hospital
%%%
%%% Manuscript
%%% 'Functional brain dynamic are shaped by connectome n-simplicial organization'
%%% Emeline Mullier, Jakub Vohryzek, Alessandra Griffa, Yasser Alemàn-Gómez, Celia Hacker, Kathryn Hess, Patric Hagmann
%%%
%%% ==========================================================================================================================



%%% ==================================================================
%%% Parameters
%%% ==================================================================
working_dir = 'DynamicMeasures';
addpath(genpath(working_dir));
load('GenerationStates/Data/CCsPrep_CCs.mat'); load('Data/Louvain_g1_init20.mat'); % load the analysis
suffixe = 'STConn_g1'; % suffixe for the name of the figures saved
fig_path = 'DynamicMeasures/Figures/'; % path to save the figure
export = 1; % save the figures (1) or not (0)
scale = 4; % scale of the parcellation
thNoise = 2;  % threshold to consider a community as noise. 
labels_clus = {'State1', 'State2', 'State3', 'State4', 'State5',  'State6', 'State7', 'State8', 'State9', 'State10', 'State11', 'State12'}; % Labels of the states (to modify according to your states)
cmap = lines; col = cmap(1:length(labels_clus),:); %% to modify according to the colors you want, here we randomly picked the colormap 'lines'
fstitle =16; fslabel =12; % fontsize for plotting
Surfl = Read_Surface('PlotKit/fsaverage/surf/lh.pial'); Surfr = Read_Surface('PlotKit/fsaverage/surf/rh.pial'); % Read the right and left hemispheres surfaces for plotting
annotr = sprintf('PlotKit/fsaverage/label/rh.lausanne2008.scale%d.annot', scale); annotl = sprintf('PlotKit/fsaverage/label/lh.lausanne2008.scale%d.annot', scale); % Path to annot files of the parcellation
load('PlotKit/labels_index_CORTICAL_Laus2008_all_scales.mat'); % load indices of the cortical regions. 


%%% ==================================================================
%%% Re-attribute the maps of noisy communities to the closest state.
%%% ==================================================================
Noise =  accumarray(CIU, ones(size(CIU)));
idxsNoise = find(Noise<thNoise);
idxsClus = [1:max(CIU)]; idxsClus(idxsNoise) = [];
idxsCIU = find(ismember(CIU,idxsNoise));
brain_maps_nonzeros = brain_maps(ixc{scale},sum(brain_maps,1) ~= 0);
mat_clusters = mat_clusters(ixc{scale}, find(Noise>thNoise));
nbStates = size(mat_clusters,2);

if export
    for i=1:nbStates
        Plot_Brain_5views(Surfr, Surfl, annotr, annotl, mat_clusters(:,i), 'colMap', 'jet');
        Figurename = [fig_path sprintf('Plot_5views_State%d_%s.tiff', i, suffixe)];
        export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r200' );
        close;
    end
end

for i=1:length(idxsCIU)
    ls_corr = [];
    for j=1:size(mat_clusters,2)
        r = corrcoef(brain_maps_nonzeros(:,idxsCIU(i)), mat_clusters(:,j));
        ls_corr = [ls_corr r(1,2)];
    end   
    CIU(idxsCIU(i)) = idxsClus(find(ls_corr == max(ls_corr)));
end

%%% Re-number the states
nbStates = length(unique(CIU));
CIU2 = zeros(size(CIU));
for i=1:length(idxsClus)
   CIU2(find(CIU==idxsClus(i)))=i;
end
CIU = CIU2;
CIU_silent = zeros(1,size(brain_maps,2));
idxsCIU = find(sum(brain_maps,1) ~= 0);
CIU_silent(idxsCIU) = CIU;
CIU_silent = CIU_silent + 1;


%%% ==================================================================
%%% Dynamic measures
%%% ==================================================================

%%% Occurence
%%%%%%%%%%%%%%%%
for st = 1:length(unique(CIU))
    pr_occurence(st) = (sum(CIU == st)./length(CIU))*100;
end
figure;
for i=1:length(pr_occurence)
    hold on;
    b = bar(i, pr_occurence(i));
    b.FaceColor = col(i,:);
end
set(gca, 'XTick', [1:max(CIU)],'XTickLabels', labels_clus(1:max(CIU)), 'fontsize',14, 'TickLabelInterpreter', 'latex');
ylabel('Occurence (Percentage)','interpreter', 'latex');
title('States occurrences', 'FontSize', 18, 'interpreter', 'latex');
xtickangle(45);
if export 
    save(sprintf('%s/node_size_%s.txt',fig_path, suffixe),'pr_occurence','-ascii')
    Figurename = [fig_path sprintf('StatesOccurenceLausanne_%s.tiff',suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
    close;
end


%%% Inter-states transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbStates_silent = length(unique(CIU_silent));
occ_transitions = zeros(1,nbStates_silent*nbStates_silent);
for i=1:length(CIU_silent)-1
    occ_transitions(sub2ind([nbStates_silent,nbStates_silent],CIU_silent(i), CIU_silent(i+1))) = occ_transitions(sub2ind([nbStates_silent,nbStates_silent],CIU_silent(i), CIU_silent(i+1))) + 1;
end
occ_transitions = tril(occ_transitions,-1)+triu(occ_transitions,1);
transfer_matrix = reshape(occ_transitions, [nbStates_silent nbStates_silent]);
transfer_matrix = transfer_matrix(2:end,2:end);
transfer_matrix = tril(transfer_matrix,-1)+triu(transfer_matrix,1);
for i=1:size(transfer_matrix,1)
    transfer_matrix(i,:) = transfer_matrix(i,:)./sum(transfer_matrix(i,:));
end

figure;
imagesc(transfer_matrix);
t = transfer_matrix; t = round(t,2);
t = num2cell(t); [x,y] = meshgrid(1:12); t = cellfun(@num2str, t, 'UniformOutput', false); 
title('Inter-States Transition Probabilities', 'FontSize',16,'interpreter', 'latex');
set(gca, 'XTick', [1:12], 'YTick', [1:12],'XTickLabels', labels_clus, 'YTickLabels', labels_clus, 'TickLabelInterpreter', 'latex');
xtickangle(45); colorbar;
if export
    Figurename = [fig_path sprintf('OccurenceTransitionsMatrix_%s.tiff',suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
    close;
    save(sprintf('%s/transfer_matrix_%s.txt', fig_path, suffixe), 'transfer_matrix', '-ascii');
end


%%% ==================================================================
%%% Circular Graph
%%% ==================================================================


%%% Save the files for the circular graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Save the probability of occurrence of each state as node size
fid = fopen(sprintf('%s/Data/Labels.csv', working_dir),'w');
for i=1:max(CIU)
    fwrite(fid, sprintf('%s \n', labels_clus{i}));
end
%%% Save the colors corresponding to the states
csvwrite(sprintf('%s/Data/Colors.csv',working_dir), col); 
%%% Save the inter-states transition matrix
csvwrite(sprintf('%s/Data/TM.csv',working_dir), transfer_matrix); 


%%% Plot the transition matrix using the Circular Graph representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl = 'manitou';
th = 0;
txt_STConn = importdata(sprintf('%s/Data/Labels.csv', working_dir));
Colors_STConn = csvread(sprintf('%s/Data/Colors.csv', working_dir));
TM = csvread(sprintf('%s/Data/TM.csv', working_dir));
Plot_Circular_CMatrix(TM, 'connThresh', th, 'stColors', Colors_STConn, 'stNames', txt_STConn,  'directed', 1);
if export
    Figurename = [fig_path sprintf('CircularGraph_%s.tiff',suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
    close;
end
  