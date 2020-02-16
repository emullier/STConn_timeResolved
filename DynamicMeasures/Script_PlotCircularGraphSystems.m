


fig_path = 'Figures/';

cl = 'manitou';
th = 1/12; %%
th = 1/48;

%txt_STConn = {'DA', 'LIMB', 'VA', 'VIS', 'FP-DMN',  'VA-SM', 'LIMB', 'FP', 'dorsal-SM', 'DMN', 'DMN-FP', 'ventral-SM'};
txt_STConn = {'DA', 'inf-LIMB', 'AUD', 'VIS', 'left-FP',  'VA', 'sup-LIMB', 'right-FP', 'dorsal-SM', 'post-DMN', 'ant-DMN', 'ventral-SM'};
Colors_STConn = [0 118 14; 220 248 164; 139 69 19; 120 18 134; 230 148 34; 196 58 250; 220 248 164; 230 148 34; 70 130 180; 205 62 78; 205 62 78; 70 130 180]./255;
txt_PPA =  {'VIS', 'FP', 'SM', 'VA', 'DMN', 'DA', 'LIMB', 'VA', 'DMN', 'FP-LIMB-DMN', 'SM', 'LIMB-FP'};
Colors_PPA = [120 18 134; 230 148 34; 70 130 180; 196 58 250; 205 62 78; 0 118 14; 220 248 164; 196 58 250; 205 62 78; 230 148 34; 70 130 180; 220 248 164];
ocurrPerc_STConn = load([fig_path 'node_size_STConn_12states.txt']); 
ocurrPerc_PPA = load([fig_path 'node_size_PPA_12states.txt']); 

orders = [1 0 0 0 1 1 0 1 0 1 1 0];
idxs = [10 11 5 8 2 7 6 1  4 3 12 9];

Colors = zeros(size(Colors_STConn,1),size(Colors_STConn,2),2);
Colors(:,:,1) = Colors_STConn/255;
Colors(:,:,2) = Colors_PPA/255;

clear txt
txt{1,:} = txt_STConn;
txt{2,:} = txt_PPA;

Ocurr = zeros(size(ocurrPerc_STConn,1),size(ocurrPerc_STConn,2),2);
Ocurr(:,:,1) = ocurrPerc_STConn;
Ocurr(:,:,2) = ocurrPerc_PPA;

nbStates = 12;
TM = zeros(nbStates,nbStates,2);
transfer_matrix = load([fig_path 'transfer_matrix_STConn_12states.txt']);
%transfer_matrix(transfer_matrix<th) = 0;
load('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/Randomization/Figures/avgTransMatPerm_nbPerm10000.mat');
load('/home/localadmin/Documents/STConn_Analyses/102019_BrainDecoding_CCs/Randomization/Figures/pValsTransMat_nbPerm10000.mat');

transfer_matrix(pVals>0.05) = 0;
transfer_matrix(transfer_matrix<avgTransMatPerm) = 0;

TM(:,:,1) = transfer_matrix;
transfer_matrix = load([fig_path 'transfer_matrix_PPA_12states.txt']);
transfer_matrix(transfer_matrix<th) = 0;
TM(:,:,2) = transfer_matrix;

TM = TM(idxs,idxs,1);
Cols = Colors_STConn(idxs,:);
txt = txt_STConn(idxs);
Occ = Ocurr(:,:,1)'; Occ = Occ(idxs);

txt2=[];
for i=1:length(txt)
   txt2{i} = sprintf('%s \n %g ', txt{i}, ceil(Occ(i))/100);
   txt2{i} = strrep(txt2{i}, '-','');
end

% Plot_Circular_CMatrix(TM, 'connThresh', th, 'stColors', Cols, 'stNames', txt, 'CharFile',Occ , 'directed', 1);
%Plot_Circular_CMatrix(TM, 'connThresh', th, 'stColors', Cols, 'stNames', txt, 'directed', 1,'hemiVar',1, 'CharFile',Occ );
%Plot_Circular_CMatrix_asym(TM, 'connThresh', th, 'stColors', Cols, 'stNames', txt);
%Plot_Circular_CMatrix_asym(TM, 'connThresh', th, 'stColors', Cols, 'stNames', txt, 'directed', 1,'hemiVar',1, 'CharFile',Occ );
%Plot_Circular_CMatrix_asym(TM, 'connThresh', th, 'stColors', Cols, 'stNames', txt, 'directed', 1,'hemiVar',1, 'CharFile', Occ, 'asym', [8,4]);
%Plot_Circular_CMatrix_asym(TM', 'connThresh', th, 'stColors', Cols, 'stNames', txt, 'directed', 1,'hemiVar',1, 'CharFile', Occ, 'asym', [8,4], 'invbool',1);


Plot_Circular_CMatrix_asym(TM, 'connThresh', th, 'stColors', Cols, 'stNames', txt2, 'directed', 1,'hemiVar',1,  'asym', [8,4]);
Figurename = [fig_path sprintf('CirGraphSystemsWhite_Direction1_th%d_%s.tiff',th, 'STConn')];
export_fig(Figurename,'-tiff', gca, '-nocrop','-transparent','-opengl','-r300' );
close;

Plot_Circular_CMatrix_asym(TM', 'connThresh', th, 'stColors', Cols, 'stNames', txt2, 'directed', 1,'hemiVar',1,  'asym', [8,4], 'invbool',1);
Figurename = [fig_path sprintf('CirGraphSystemsWhite_Direction2_th%d_%s.tiff',th, 'STConn')];
export_fig(Figurename,'-tiff', gca, '-nocrop','-transparent','-opengl','-r300' );
close;


%Plot_Circular_CMatrix_asym(TM', 'invbool', 1,'connThresh', th, 'stColors', Cols, 'stNames', txt, 'directed', 1,'hemiVar',1, 'CharFile', Occ);


% methods = {'STConn_12states', 'PPA_12states'}
% for i=1:length(methods)
%     Plot_Systems_Graph(TM(:,:,i), cl, Ocurr(1,:,i),0, txt{i,:}, Colors(:,:,i));
%     Figurename = [fig_path sprintf('CirGraphSystemsWhite_Direction1_th%d_%s.tiff',th,methods{i})];
%     export_fig(Figurename,'-tiff', gca, '-nocrop','-transparent','-opengl','-r300' );
%     close;
%     Plot_Systems_Graph(TM(:,:,i)', cl, Ocurr(1,:,i),1, txt{i,:}, Colors(:,:,i));
%     Figurename = [fig_path sprintf('CirGraphSystemsWhite_Direction2_th%d_%s.tiff',th, methods{i})];
%     export_fig(Figurename,'-tiff', gca, '-nocrop','-transparent','-opengl','-r300' );
%     close;
% end


