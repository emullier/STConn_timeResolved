



function varargout = Plot_Brain_5views(Surfr, Surfl, annotr, annotl, values, varargin);

% Syntax :
%  Plot_Brain_5views(Surfr, Surfl, annotr, annotl, values, varargin);
%
% This function plots regionwise values on a brain with 5 views:
% lateral left and right, sagittal left and right and superior. 
%
% Input Parameters:
%   Surfr, Surfl       : Surfaces files of right anf left hemispheres
%   annotr, annotl     : Annot files of right and left hemispheres of the
%                        chosen parcellation (used for the regionwise values)
%   values             : List of the regionwise values
%
% Optional Parameters:
%   transpVal          : Transparency vector(values have to be between 0 and 1).
%   colMap             : Colormap used to see the results.
%   figcolor           : Figure color  
%   inflated           : 0 if pial surface, 1 if inflated (freesurfer surfaces)
%
%_______________________________________
% Authors: Emeline Mullier, Jakub Vohryzek, Yasser Aleman Gomez
% Lausanne University Hospital
% October 4th, 2018
% Version $1.0



    %% Input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Default values for optional parameters
    if nargin<5 
        error('5 inputs are mandatory: Surfr, Surfl, annotr, annotl, values');
    else
        transpVal = 1;      % Opacity
        colMap = 'winter';  % ColorMap
        newFig = 'y';       % New figure
        % Figure Properties
        param.figcolor = [1 1 1];     % Color
        param.figurevisible = 'on';   % Visible
        param.newfigure = 1    ;      % Boolean variable to create new figure
        set(0,'units','centimeters');
        cm_screen = get(0,'screensize');
        figPosition = [1 1 cm_screen(3)-3 cm_screen(4)-3*cm_screen(3)/cm_screen(4)]; % Position
        figUnits = 'centimeters';   
        inflated = 0; % Units
    end

    % If optional arguments are provided
    if numel(varargin)>1
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'transpVal' % Group Separation
                    transpVal=varargin{2};
                case 'colMap' % Jitter compression percent
                    colMap=varargin{2};
                case 'figcolor'
                    param.figcolor=varargin{2};
                case 'figurevisible'
                    param.figurevisible=varargin{2};
                case 'figUnits'
                    figUnits=varargin{2};
                case 'figPosition'
                    figPosition=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                case 'inflated'
                    inflated=varargin{2};   
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
                end
                varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end

    
    %% Attribution of the region values to the vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Yeo's 7FS
    %colMap = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]/250;
    %%% Colors Lobes
    %hColor = jet; 
    %colMap = hColor(1:7:63,:);
    %labelFS = {'unknown','Frontal','Corpus', 'Cingulate', 'Occipital', 'Temporal', 'Parietal', 'Insula', 'Subcortical'};    
   

 
    % Read annot files for our parcellation
    [txtr,ctabr,colorsr] = read_cfiles(annotr);
    [txtl,ctabl,colorsl] = read_cfiles(annotl);

    % Apply the parcellation to the vertices
    Surfr.SurfData.FaceVertexCData = colorsr;
    Surfl.SurfData.FaceVertexCData = colorsl;
    
    % Match the vertices and the regionwise values - right hemisphere
    ind_ctxr = [];
    Colors = zeros(length(txtr),3);
    Is = zeros(length(txtr),1);
    for i = 1:size(ctabr.table,1)-1
        ind = find(txtr == ctabr.table(i+1,5));
        Is(ind) = values(i); 
        %Colors(ind,:) = repmat(colMap(values(i),:), size(ind));
        ind_ctxr = [ind_ctxr; ind];   
    end
    % Rescale the values on the colormap
    Surfr.Is = Is; Colors = Surf_Color(Surfr,colMap); Colors(Is==0,:) = 1 ;
    Surfr.SurfData.FaceVertexCData = Colors; 
    % Changes subcortical vertices to black
    ind_subctxr = find(~ismember([1:length(txtr)], ind_ctxr));
    Surfr.SurfData.FaceVertexCData(ind_subctxr,:) = zeros(size(Surfr.SurfData.FaceVertexCData(ind_subctxr,:))); 

        
    % Match the vertices and the regionwise values - left hemisphere
    ind_ctxl = [];
    Colors = zeros(length(txtl),3);
    Is = zeros(length(txtl),1);
    for i = 1:size(ctabl.table,1)-1
        ind = find(txtl == ctabl.table(i+1,5));
        Is(ind) = values(i+size(ctabr.table,1)-1);
        ind_ctxl = [ind_ctxl; ind];   
        %Colors(ind,:) = repmat(colMap(values(i),:), size(ind));
    end
    Surfl.Is = Is; Colors = Surf_Color(Surfl,colMap); 
    Colors(Is==0,:) = 1 ; Surfl.SurfData.FaceVertexCData = Colors; 
    ind_subctxl = find(~ismember([1:length(txtl)], ind_ctxl));
    Surfl.SurfData.FaceVertexCData(ind_subctxl,:) = zeros(size(Surfl.SurfData.FaceVertexCData(ind_subctxl,:))); 
    if inflated
        Surfl.SurfData.vertices(:,1,1) = Surfl.SurfData.vertices(:,1,1)-65;
    end

    % Merge right and left hemispheres surfaces
    Surf = Compound_Surf([Surfr; Surfl]);
    
    
    %% Plot the 5 views brain visualization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create the structure for the 5 views visualization
    Surf5{1,1} = Surfr; Surf5{2,1} = Surfr; Surf5{3,1} = Surf; Surf5{4,1} = Surfl; Surf5{5,1} = Surfl;
    Nsurf = length(Surf5);
    list_views = {[90 0], [-90 0], [0 90], [90 0], [-90 0], [90 90 -180],[90 0], [-90 0], [0 90], [90 0], [-90 0], [90 90 -180]};
    
    FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
    
    % Subplot of each surface
    for i = 1:Nsurf 
    
        temp = Surf5{i,1};
        clmap = colMap;
        axIds(i) = subplot_tight(1, Nsurf,i,[0.03]);
        %axIds(i) = subplot_tight(Nsurf,1,i,[0.05]);
                
        if ischar(temp)
            if exist(deblank(temp(1,:)),'file')
                for j = 1:size(temp,1);
                    tempSurf = deblank(temp(j,:));
                    Surf(j) = Surface_Checking(tempSurf);
                end
                [hSurf, jSurf, surfRealMap] = test_Plot_Surf(temp,'FigID',FigID,'transpVal',1,'colMap',clmap);grid on; 
            else
                error('Wrong surface filename');
                return
            end
        else
            [hSurf, jSurf, surfRealMap] = test_Plot_Surf(temp,'FigID',FigID,'transpVal',1,'colMap',clmap);grid on;
        end
        
        view(axIds(i),list_views{i});
        axis off;
        colorbar('off');
        
        % Plot the colorbar on the last view
        if i == Nsurf
            colormap(colMap);
            H = colorbar('Southoutside');
            H.Position(2) = H.Position(2) + 0.25
        end
         
        
    end

end