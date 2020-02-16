

cd '/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/DynamicMeasures';

fig_path = 'Figures/';

cl = 'manitou';
%th = 1/12; %%
th = 1/24;

HCP = 2;

labels_yeo = {'VIS', 'SM', 'DA', 'VA', 'LIMB', 'FP', 'DMN'};
col_yeo = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./255;
order_yeo = [0 0 1 1 0 1 1];
col = []; 
if HCP == 1
    name = '50HCP1';
    labels_clus = {'VIS-DA', 'SM', 'VA', 'DMN-FP1', 'DA', 'DMN-FP2', 'vSM', 'FP', 'dSM', 'DMN', 'VIS'};
    idx_yeo = [1 2 4 7 3 7 2 6 2 7 1];
    for i=1:length(idx_yeo)
        id = idx_yeo(i);
        col = [col; col_yeo(id,:)]; 
    end
    Ocurr = load(sprintf('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/DynamicMeasures/Figures/node_size_%s_12states.txt',name));
    nbStates = length(labels_clus);
    TM = zeros(nbStates,nbStates,2);
    transfer_matrix = load([fig_path 'transfer_matrix_' name '_12states.txt']);
    transfer_matrix(transfer_matrix<th) = 0;
elseif HCP == 2
    name = '50HCP2';
    labels_clus = {'DMN-FP', 'DMN', 'VA1', 'VIS', 'FP', 'v-SM1', 'VA2', 'FP-VA', 'SM-VIS', 'v-SM2', 'VIS-DA', 'dSM'};
    idx_yeo = [7 7 4 1 6 2 4 6 2 2 1 2];
    for i=1:length(idx_yeo)
        id = idx_yeo(i);
        col = [col; col_yeo(id,:)]; 
    end
    Ocurr = load(sprintf('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/DynamicMeasures/Figures/node_size_%s_12states.txt',name));
    nbStates = length(labels_clus);
    TM = zeros(nbStates,nbStates,2);
    transfer_matrix = load([fig_path 'transfer_matrix_' name '_12states.txt']);
    transfer_matrix(transfer_matrix<th) = 0;
end


Plot_Systems_Graph(transfer_matrix, cl, Ocurr,0, labels_clus, col);
Figurename = [fig_path sprintf('CirGraphSystemsWhite_Direction1_th%d_%s.tiff',th,name)];
export_fig(Figurename,'-tiff', gca, '-nocrop','-transparent','-opengl','-r300' );
close;
Plot_Systems_Graph(transfer_matrix', cl, Ocurr,1, labels_clus, col);
Figurename = [fig_path sprintf('CirGraphSystemsWhite_Direction2_th%d_%s.tiff',th, name)];
export_fig(Figurename,'-tiff', gca, '-nocrop','-transparent','-opengl','-r300' );
close;



