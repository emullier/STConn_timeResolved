

%%%
%%%
%%% Updated version without integrating the disoverlaps.
%%% An extra step is required to attribute the time points attributed to
%%% noisy community to a state 
%%%
%%% EM
%%% 07.10.2019
%%% Lausanne University Hospital
%%%


%addpath(genpath('/home/localadmin/Documents/STConn_Analyses'));

cd '/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/GenerationStates';
working_dir = '/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/DynamicMeasures';

load('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/GenerationStates/Data/CCsPrep_CCs.mat');
load('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/GenerationStates/Data/Louvain_g2_init20.mat');
suffixe = 'STConn_12states';

%load('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/GenerationStates/Data/CCsPrep_CCs_50HCP.mat');
%load('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/GenerationStates/Data/Louvain_g2_init20_50HCP.mat');
%suffixe = '50HCP_11states';

cd '/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/DynamicMeasures';
fig_path = 'Figures/';

scale = 4; thNoise = 50; fstitle =16; fslabel =12;
export = 0;

%%% TO BE CHECKED
labels_clus = {'DA', 'inf-LIMB', 'AUD', 'VIS', 'left-FP',  'VA', 'sup-LIMB', 'right-FP', 'dorsal-SM', 'post-DMN', 'ant-DMN', 'ventral-SM'};
col = [0 118 14; 220 248 164; 139 69 19; 120 18 134; 230 148 34; 196 58 250; 220 248 164; 230 148 34; 70 130 180; 205 62 78; 205 62 78; 70 130 180]./255;
order_system = [1 1 0 0 1 1 1 1 0 1 1 0]; 


% Read the right and left hemispheres surfaces
Surfl = Read_Surface('/opt/freesurfer/subjects/fsaverage/surf/lh.pial');
Surfr = Read_Surface('/opt/freesurfer/subjects/fsaverage/surf/rh.pial');
% Path to annot files of the parcellation
annotr = sprintf('/opt/freesurfer/subjects/fsaverage/label/rh.lausanne2008.scale%d.annot', scale);
annotl = sprintf('/opt/freesurfer/subjects/fsaverage/label/lh.lausanne2008.scale%d.annot', scale);
load('/mnt/data/RAD/HAGMANN_GROUP/EMullier/MATLAB/EEG_PROJECT/plot/plot_kit/needed_data/labels_index_CORTICAL_Laus2008_all_scales.mat');


%%% ==================================================================
%%% Re-attribute the maps of noisy communities to the closest state.
%%% ==================================================================
Noise =  accumarray(CIU, ones(size(CIU)));
idxsNoise = find(Noise<thNoise);
idxsClus = [1:max(CIU)]; idxsClus(idxsNoise) = [];
idxsCIU = find(ismember(CIU,idxsNoise));
brain_maps_nonzeros = brain_maps(ixc{scale},sum(brain_maps,1) ~= 0);
mat_clusters = mat_clusters(ixc{scale}, find(Noise>thNoise));
%brain_maps_nonzeros = brain_maps(:,sum(brain_maps,1) ~= 0); %% for the case HCP
%mat_clusters = mat_clusters(:, find(Noise>thNoise)); %% for the case HCP
nbStates = size(mat_clusters,2);

% if export
%     for i=1:nbStates
%         Plot_Brain_5views(Surfr, Surfl, annotr, annotl, mat_clusters(:,i), 'colMap', 'jet');
%         Figurename = [fig_path sprintf('Plot_5views_State%d_%s.tiff', i, suffixe)];
%         export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r200' );
%         close;
%     end
% end



%%% Solution1: Include a noisy state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIU(idxsCIU) = 0; %%% Put the time points attributed to noisy communities to zero
% 
% %%% Re-number the states in case of one of the noisy clus was among the first ones. 
% CIU2 = zeros(size(CIU));
% for i=1:length(idxsClus)
%     CIU2(find(CIU==idxsClus(i)))=i;
% end
% CIU = CIU2;
% 
% %%% change the numbering to have the noisy time points set to 13 (nbStates+1) instead of 0
% CIU(find(CIU==0)) = nbStates + 1; %% 13 = noisy states
% 
% %%% include the silent time points with 0 state
% CIU_silent = zeros(1,size(brain_maps,2));
% idxsCIU = find(sum(brain_maps,1) ~= 0);
% CIU_silent(idxsCIU) = CIU;
% CIU_silent = CIU_silent + 1;


%%% Solution 2: Attribution of the clusters according to their closest spatially
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(idxsCIU)
    ls_corr = [];
    for j=1:size(mat_clusters,2)
        r = corrcoef(brain_maps_nonzeros(:,idxsCIU(i)), mat_clusters(:,j));
        ls_corr = [ls_corr r(1,2)];
    end   
    CIU(idxsCIU(i)) = idxsClus(find(ls_corr == max(ls_corr)));
end

%%% Re-number the states in case of one of the noisy clus was among the first ones. 
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


%%% ================ %%%
%%% DYNAMIC ANALYSIS %%%
%%% ================ %%%

idxs = [10 11 5 8 2 7 6 1  4 3 12 9];
order_system = order_system(idxs)
%%% Occurence
for st = 1:length(unique(CIU))
    pr_occurence(st) = (sum(CIU == st)./length(CIU))*100;
end
pr_occurence = pr_occurence(idxs);
col = col(idxs,:);
labels_clus = labels_clus(idxs);
figure;
for i=1:length(pr_occurence)
    hold on;
    b = bar(i, pr_occurence(i));
   b.FaceColor = col(i,:);
end
set(gca, 'XTick', [1:12],'XTickLabels', labels_clus, 'fontsize',14, 'TickLabelInterpreter', 'latex');
ylabel('Occurence (Percentage)','interpreter', 'latex');
title('States occurrences', 'FontSize', 18, 'interpreter', 'latex');
xtickangle(45);
%csvwrite(sprintf('%s/%snode_size_%s.csv',working_dir, fig_path, suffixe),pr_occurence)
if export 
    save(sprintf('%s/%snode_size_%s.txt',working_dir, fig_path, suffixe),'pr_occurence','-ascii')
    Figurename = [fig_path sprintf('StatesOccurenceLausanne_%s.tiff',suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
    close;
end

nbStates_silent = length(unique(CIU_silent));
occ_transitions = zeros(1,nbStates_silent*nbStates_silent);
for i=1:length(CIU_silent)-1
    occ_transitions(sub2ind([nbStates_silent,nbStates_silent],CIU_silent(i), CIU_silent(i+1))) = occ_transitions(sub2ind([nbStates_silent,nbStates_silent],CIU_silent(i), CIU_silent(i+1))) + 1;
end
occ_transitions = tril(occ_transitions,-1)+triu(occ_transitions,1);
transfer_matrix = reshape(occ_transitions, [nbStates_silent nbStates_silent]);
transfer_matrix = transfer_matrix(2:end,2:end);
%transfer_matrix = transfer_matrix./sum(transfer_matrix(:));
transfer_matrix = tril(transfer_matrix,-1)+triu(transfer_matrix,1);
for i=1:size(transfer_matrix,1)
    transfer_matrix(i,:) = transfer_matrix(i,:)./sum(transfer_matrix(i,:));
end
transfer_matrix = transfer_matrix(idxs,idxs);

figure;
imagesc(transfer_matrix);
t = transfer_matrix; t = round(t,2);
t = num2cell(t); [x,y] = meshgrid(1:12); t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
%text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'interpreter', 'latex');
title('Inter-States Transition Probabilities', 'FontSize',16,'interpreter', 'latex');
set(gca, 'XTick', [1:12], 'YTick', [1:12],'XTickLabels', labels_clus, 'YTickLabels', labels_clus, 'TickLabelInterpreter', 'latex');
xtickangle(45); colorbar;
%csvwrite(sprintf('%s/%stransfer_matrix_%s.csv', working_dir, fig_path, suffixe), transfer_matrix);
if export
    Figurename = [fig_path sprintf('OccurenceTransitionsMatrix_%s_withDiag2.tiff',suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
    close;
    save(sprintf('%s/%stransfer_matrix_%s.txt', working_dir, fig_path, suffixe), 'transfer_matrix', '-ascii');
end

low = find(order_system==0);
high = find(order_system==1);

figure; 
subplot(2,2,1); imagesc(transfer_matrix(low,low)); axis square; title('Low - Low', 'fontsize',fstitle,'interpreter', 'latex'); colorbar;
set(gca, 'XTick',[1:length(low)],'XTickLabel', labels_clus(low), 'YTick',[1:length(low)],'YTickLabel', labels_clus(low), 'TickLabelInterpreter', 'latex');
subplot(2,2,2); imagesc(transfer_matrix(low,high)); axis square; title('Low - High', 'fontsize',fstitle,'interpreter', 'latex'); colorbar;
set(gca, 'XTick',[1:length(high)],'XTickLabel', labels_clus(high), 'YTick',[1:length(low)],'YTickLabel', labels_clus(low), 'TickLabelInterpreter', 'latex');
subplot(2,2,3); imagesc(transfer_matrix(high,low)); axis square; title('High - Low', 'fontsize',fstitle,'interpreter', 'latex'); colorbar;
set(gca, 'XTick',[1:length(low)],'XTickLabel', labels_clus(low), 'YTick',[1:length(high)],'YTickLabel', labels_clus(high), 'TickLabelInterpreter', 'latex');
subplot(2,2,4); imagesc(transfer_matrix(high,high)); axis square; title('High - High', 'fontsize',fstitle,'interpreter', 'latex'); colorbar;
set(gca, 'XTick',[1:length(high)],'XTickLabel', labels_clus(high), 'YTick',[1:length(high)],'YTickLabel', labels_clus(high), 'TickLabelInterpreter', 'latex');
if export
    Figurename = [fig_path sprintf('TransitionsMatrix_LowHighOrders_%s.tiff',suffixe)];
    export_fig(Figurename,'-tiff', gcf, '-nocrop','-transparent','-opengl','-r300' );
    close;
end



% TR = 1.92;
% % SWITCHING FREQUENCY
% switchFreq = mean(diff(CIU)~=0)/TR; % same measure for all clusters
% 
% clear clusProb; clear lifeTime;
% for cl = 1:length(unique(CIU))
%     % 3. CLUSTER PROBABILITY
%     clustProb(cl) = mean(CIU == cl); % Like this already normalised!
%     % 4. CLUSTER LIFETIME
%     clustDur = [];
%     clustLifeTime = [];
%     clustDur = (CIU == cl);
% 
%     % Detect switches in and out of this state
%     onset = find(diff([0; clustDur; 0]') == 1);
%     offset = find(diff([0; clustDur; 0]') == -1);
% 
%     if ~isempty(onset) && ~isempty(offset)
%         clustLifeTime = offset - onset;
%     else
%         clustLifeTime = 0;
%     end
% 
%     lifeTime(cl) = mean(clustLifeTime);
% end
% figure; bar(lifeTime);
% ylabel('Average life time (tps)');
% 




  