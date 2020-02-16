
function varargout = Plot_Brains_1view(Surfr, Surfl, annotr, annotl, matrix_values, varargin);

% Syntax :
%  Plot_Brain_5views(Surfr, Surfl, annotr, annotl, matrix_values, varargin);
%
% This function plots different regionwise values on brains  with 1 specific view to compare them. 
% For clusters comparison for example
%
% Input Parameters:
%   Surfr, Surfl       : Surfaces files of right anf left hemispheres
%   annotr, annotl     : Annot files of right and left hemispheres of the
%                        chosen parcellation (used for the regionwise values)
%   matrix_values      : NxM matrix with N number of subjects, list of
%                        values and M number of regions of the parcellation
%
% Optional Parameters:
%   transpVal          : Transparency vector(values have to be between 0 and 1).
%   colMap             : Colormap used to see the results.
%   figcolor           : Figure color  
%   view_orient        : Orientation for the brain view
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
        colMap = 'autumn';  % ColorMap
        newFig = 'y';       % New figure
        view_orient = [0 90];
        % Figure Properties
        param.figcolor = [1 1 1];     % Color
        param.figurevisible = 'on';   % Visible
        param.newfigure = 1    ;      % Boolean variable to create new figure
        set(0,'units','centimeters');
        cm_screen = get(0,'screensize');
        figPosition = [1 1 cm_screen(3)-3 cm_screen(4)-3*cm_screen(3)/cm_screen(4)]; % Position
        figUnits = 'centimeters';                                                    % Units
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
                case 'view_orient'
                    view_orient = varargin{2};
                case 'figurevisible'
                    param.figurevisible=varargin{2};
                case 'figUnits'
                    figUnits=varargin{2};
                case 'figPosition'
                    figPosition=varargin{2};
                case 'FigID'
                    FigID=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
                end
                varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end

    
    %% Attribution of the different region values to the vertices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Read annot files for our parcellation
    [txtr,ctabr,colorsr] = read_cfiles(annotr);
    [txtl,ctabl,colorsl] = read_cfiles(annotl);

    % Apply the parcellation to the vertices
    Surfr.SurfData.FaceVertexCData = colorsr;
    Surfl.SurfData.FaceVertexCData = colorsl;
    
    
    for t=1:size(matrix_values,1)
    
        % Match the vertices and the regionwise values - right hemisphere
        ind_ctxr = [];
        Is = zeros(length(txtr),1);
        for i = 1:size(ctabr.table,1)-1
            ind = find(txtr == ctabr.table(i+1,5));
            Is(ind) = matrix_values(t,i); 
            ind_ctxr = [ind_ctxr; ind];   
        end
        % Rescale the values on the colormap
        Surfr.Is = Is; Colors = Surf_Color(Surfr,colMap); Colors(Is==0,:) = 1 ;Surfr.SurfData.FaceVertexCData = Colors; 
        % Changes subcortical vertices to black
        ind_subctxr = find(~ismember([1:length(txtr)], ind_ctxr));
        Surfr.SurfData.FaceVertexCData(ind_subctxr,:) = zeros(size(Surfr.SurfData.FaceVertexCData(ind_subctxr,:))); 

        % Match the vertices and the regionwise values - left hemisphere
        ind_ctxl = [];
        Is = zeros(length(txtl),1);
        for i = 1:size(ctabl.table,1)-1
            ind = find(txtl == ctabl.table(i+1,5));
            Is(ind) = matrix_values(t,i+size(ctabr.table,1)-1);
            ind_ctxl = [ind_ctxl; ind];   
        end
        Surfl.Is = Is; Colors = Surf_Color(Surfl,colMap);Colors(Is==0,:) = 1 ; Surfl.SurfData.FaceVertexCData = Colors; 
        ind_subctxl = find(~ismember([1:length(txtl)], ind_ctxl));
        Surfl.SurfData.FaceVertexCData(ind_subctxl,:) = zeros(size(Surfl.SurfData.FaceVertexCData(ind_subctxl,:))); 

        % Merge right and left hemispheres surfaces
        Surf = Compound_Surf([Surfr; Surfl]);
        
        Surfa{t,1} = Surf;
    
    end
    
    
    %% Plot the different brains
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nsurf = length(Surfa);
    
    FigID = figure('numbertitle','off','name','New Figure','Color',param.figcolor, 'units',figUnits,'Position',figPosition,'Visible',param.figurevisible,'InvertHardcopy','off');
    
    if size(matrix_values,1) == 1
        subplotFactor = 1;
    else
        subplotFactor = 1.2 - (1/(size(matrix_values,1)));
    end

    
    % Subplot of each surface
    for i = 1:Nsurf 
    
        temp = Surfa{i,1};
        clmap = colMap;
        %axIds(i) = subplot_tight(1,Nsurf,i,[0.03]);
        axIds(i) = subplot(1,Nsurf,i);
        
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
        
        view(axIds(i),view_orient);
        axis off;
        colorbar('off');
        
        % Plot the colorbar on the last view
%         if i == Nsurf
%             H = colorbar('Southoutside');
%             H.Position(2) = H.Position(2) + 0.25
%         end
        
        %axIds(i).Position(3:4) =  hsagaxis.Position(3:4).*subplotFactor;
                  
    end
    
    
end