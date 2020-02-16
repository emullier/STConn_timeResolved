function Plot_Circular_CMatrix_asym(CMatrix, varargin);
%
% Syntax :
% Plot_Circular_Graph(CMatrix);
%
% This script plots a "bullseye" according to a specified connectivity
% matrix. Different connection colors means different networks.
%
% Input Parameters:
%       CMatrix                 : Connectivity Matrix. It can be a textfile
%       CharFile                : Cell array containing
%       param                   : Plotting parameters
%                                param.linethick - Lines thickness
%                                param.connecthr - Connectivity threshold
%
%                                Rings parameters
%                                param.N - Patches resolution
%                                param.rho - Interior and exterior rings radius
%                                param.tht - Rings theta limits
%                                param.tht0  =0;
%                                param.label =0;
%                                param.lines =[1];
%                                param.ringwidth - Ring Width
%                                param.intringwidth - Inter ring distance
%                                param.blacklines - Boolean variable to specified colored or black lines
%                                param.figcolor - Figure color
%                                param.textspace - Text space for multiple colorbars
%
% Output Parameters:
%
%
% See also: Examples proposed at line 40.
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2013
% Version $1.0

% CMatrix = '/media/PENDRIVE/scripts/Pvalue_Matrix-aparc+aseg.txt';
% CharFile = '/media/PENDRIVE/scripts/Thick+Area.txt';

% --------------------- Some Examples ------------------------------------
temp = which('Plot_Circular_Graph.m');

% Textfile - Connectivity Matrix and Characteristics File
% [pth, nm, ext] = fileparts(temp);                         % Uncomment this line to see this example
% CharFile = [pth filesep 'Thick+Area.txt'];                % Uncomment this line to see this example
% CMatrix = [pth filesep 'Pvalue_Matrix-aparc+aseg.txt'];   % Uncomment this line to see this example

% -------------------- Just Connectivity Matrices -------------------------
% ----- No 1 Single connectivity matrix
% CMatrix = zeros(8,8);      % Uncomment this line to see this example
% CMatrix(1,[2 3]) = 0.5;    % Uncomment this line to see this example
% CMatrix(4,[5 8]) = 1;      % Uncomment this line to see this example
% CMatrix([2 4 5 6],8) = 1;  % Uncomment this line to see this example
% C = (CMatrix+CMatrix')/2;  % Uncomment this line to see this example

% ----- No 2 Single connectivity matrix but with a predefined structure organization
% CMatrix = zeros(68,68);  % Uncomment this line to see this example
% CMatrix(33,56) = 1;      % Uncomment this line to see this example
% CMatrix(28,54) = 1;      % Uncomment this line to see this example
% CMatrix(54,28) = 1;      % Uncomment this line to see this example
% CMatrix(26,54) = 1;      % Uncomment this line to see this example
% CMatrix(54,26) = 1;      % Uncomment this line to see this example
% CMatrix(34,68) = 1;      % Uncomment this line to see this example
% CMatrix(68,34) = 1;      % Uncomment this line to see this example

% Textfile - Connectivity Matrix Only
% CMatrix = [pth filesep '0007-20091002-Connectivity_Matrix_FAwACP-aparc+aseg.txt'];  % Uncomment this line to see this example



%% ====================== Checking input parameters ===================== %
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    % Default Structure Names
    stNames = cellstr(num2str((1:length(CMatrix))','%03d'));
    
    % Default connectivity value
    connThresh = 0.7;             % Connectivity threshold
    
    % Default edge style
    edgeStyle = 'network'; % It will take the default edgeColors
    
    % Default edge colors
    edgeColors = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    
    % If matrix directed, gradient of colors for the edges and colormap
    directed = 0;
    invbool = 0;
    hemiVar = 0;
    idxs = [8,4];
    
    
    % Default Region colors
    stColors = [ 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9;0 1 0;];
    %stColors1 = [240,163,255;0,117,220;153,63,0;76,0,92;0,92,49;43,206,72;...
    %    255,204,153;148,255,181;143,124,0;157,204,0;194,0,136;...
    %    0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
    %    224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5]/255;
    %stColors = [stColors;stColors1];
    
    % Plotting parameters
    param.linethick = [0.4 5];      % Lines thickness
    param.blacklines = 0;         % Boolean variable to specified colored or black lines
    param.figcolor = 'white';     % Figure color
    param.textspace = 0.01;       % Text space for multiple colorbars
    param.newfigure = 1    ;       % Boolean variable to create new figure
    
    % Rings parameters
    param.N     = 10;             % Patches resolution
    param.rho   =[1 1.1];           % Interior and exterior rings radius
    param.tht   =[0 360];         % Rings theta limits
    param.tht0  =0;
    param.label =0;
    param.lines =[1];
    param.ringwidth = 0.07;       % Ring Width
    param.intringwidth = 0.05;   % Inter ring distance
end

% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('One input is mandatory');
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}
                case 'connThresh' % Connectivity Threshold
                    connThresh=varargin{2};
                case 'edgeStyle' % Edge color style (network, black or a colormap (ie. hot, jet))
                    edgeStyle=varargin{2};
                case 'CharFile' % Characteristic Map
                    CharFile=varargin{2};
                case 'stNames' % Structures Names
                    stNames =varargin{2};
                case 'param' % Default parameters for plotting
                    param=varargin{2};
                case 'stColors' % Colors for the GM structures
                    stColors =varargin{2};
                case 'edgeColors' % Colors for Edges
                    edgeColors =varargin{2};
                case 'hemiVar' % Boolean variable to separate both hemispheres
                    hemiVar =varargin{2};
                case 'directed' % gradient of colors for the edges
                    directed = varargin{2};
                case 'invbool' % in the directed case, boolean to change the direction of the gradient (Be careful: the CMatrix needs to be transposed in this case)
                    invbool =varargin{2};
                case 'asym'
                    idxs = varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

if ~isfield(param,'linethick')
    param.linethick = [0.1 5];      % Lines thickness
end
if ~isfield(param,'blacklines')
    param.blacklines = 0;         % Boolean variable to specified colored or black lines
end
if ~isfield(param,'figcolor')
    param.figcolor = 'white';     % Figure color
end
if ~isfield(param,'textspace')
    param.textspace = 0.01;       % Text space for multiple colorbars
end
if ~isfield(param,'newfigure')
    param.newfigure = 0   ;       % Boolean variable to create new figure
end

% Rings parameters
if ~isfield(param,'N ')
    param.N     = 10;             % Patches resolution
end
if ~isfield(param,'rho')
    param.rho   =[1 1.1];           % Interior and exterior rings radius
end
if ~isfield(param,'tht')
    param.tht   =[0 360];         % Rings theta limits
end
if ~isfield(param,'tht0')
    param.tht0  =0;
end
if ~isfield(param,'label')
    param.label =0;
end
if ~isfield(param,'lines')
    param.lines =[1];
end
if ~isfield(param,'ringwidth')
    param.ringwidth = 0.15;       % Ring Width
end
if ~isfield(param,'intringwidth')
    param.intringwidth = 0.025;   % Inter ring distance
end

%% ===================== Checking input parameters ====================== %
% close all;


%% ================== End of Checking input parameters ================== %

%% ================= Colormaps for multiple colorbars =================== %
cl_left = 'blackred';
cl_right= 'blackgreen';
cl_inter = 'blackblue';

colorstrings = strvcat('winter(256)','bone(256)','hsv(256)','jet(256)','pink(256)','winter(256)','gray(256)');
colorstrings = strvcat('colormaps_colors(''saginaw'',256)',...
    'colormaps_colors(''redgreen'',256)',...
    'colormaps_colors(''munising'',256)',...
    'colormaps_colors(''bluered'',256)',...
    'colormaps_colors(''pink'',256)',...
    'colormaps_colors(''jet'',256)');

%% ================ End of Colormaps for multiple colorbars ============= %

% ---------  XY coordinates of nodes in a circular layout  ---------------%


N = size(CMatrix,1);
ang = linspace(0+eps,2*pi-eps,N+1); ang = ang(2);
%ang = CharFile/(2*pi+eps); 
switch hemiVar
    case 1  
        %theta = zeros(N,1);
        %theta(1) = pi/2 + ang(1);
        %for i=2:N
        %    theta(i) = theta(i-1) + ang(i);
        %end
        theta = linspace(pi/2+ang,3/2*pi-ang,(idxs(1)))'; %theta(end) = [];
        theta2 = linspace(pi/2-ang,-pi/2+ang,(idxs(2)))'; %theta2(end) = [];   
        theta = [theta;theta2];
    otherwise
        theta = linspace(pi/2,5/2*pi,N+1)'; theta(end) = [];
end

xy = param.rho(1).* [cos(theta) sin(theta)];

% -------------------- Network Labelling -------------------------------- %

% -------------------- Network Colors ----------------------------------- %
%CMatrix = (CMatrix+CMatrix')/2;
LabNet = Label_Graph_Components(CMatrix, connThresh);
LabNet = triu(LabNet);
NNetworks = unique(nonzeros(LabNet(:)));
NN = length(NNetworks);
Ncolor = size(edgeColors,1);

re = floor(NN/Ncolor);
Col = repmat(edgeColors,[re+1 1]);
edgeColors = Col(LabNet(find(LabNet)),:);
if param.blacklines
    switch param.figcolor
        case 'black'
            edgeColors = edgeColors*0+1;
        case 'white'
            edgeColors = edgeColors*0;
    end
end


% -------------------- Points organization ------------------------------ %
ind = find(LabNet ~=0);
[xrow,xcol] = ind2sub(size(LabNet),ind);
nxrow = xrow;
nxcol = xcol;
inddel = find((xcol-xrow)>N/2); % Detecting initial point
nxrow(inddel) = xcol(inddel);
nxcol(inddel) = xrow(inddel);
xrow = nxrow;
xcol = nxcol;

% -------------------- Lines Thickness ---------------------------------- %
[ut,i1t,it2] = unique(CMatrix(ind));
LG = linspace(param.linethick(1),param.linethick(2),length(ut))';
LG = LG(it2);

% ------------------- Edge Colors
switch edgeStyle
    case 'black'
        edgeColors = edgeColors*0+1;
    case 'network'
    otherwise
        edgeColors = Val2colors(CMatrix(ind),edgeStyle);
end


% -----------------------  Bisector radius  ------------------------------%
b = round(dist(xy')*10000)/10000;
[ut,i1t,it2] = unique(b(:));
bis_rad = linspace(0.8*param.rho(1),0.1*param.rho(1),length(ut));
r_bis = reshape(bis_rad(it2),[size(b)]);
r_bis = r_bis(ind);

% -----------------------  Bisector angles  ------------------------------%
theta_bis = theta(xrow) - (theta(xrow)-theta(xcol))/2;

% ------------------- XY Bisector Coordinates  ---------------------------%
xy_bis = [r_bis r_bis].* [cos(theta_bis) sin(theta_bis)];

if param.newfigure
    switch param.figcolor
        case 'black'
            colordef black; hf = figure('numbertitle','off','name','Circular Graph','Color',[0 0 0],'Position',[0 0 1200 900]);
        case 'white'
            colordef white; hf = figure('numbertitle','off','name','Circular Graph','Color',[1 1 1],'Position',[0 0 1200 900]);
    end
end
axis off;axis equal;
Nstruct = size(LabNet,1);
if sum(LabNet(:))~=0
    % ----------------- Points to create the parabolas -----------------------%
    Par2plot = [xy(xrow,:) xy_bis xy(xcol,:) ];
    hold on
    MaxLG = max(LG);
    %Colors = repmat([1,1,1],20,1);
    
    % t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
    [row,col, v] = find(LabNet);
    for i = 1:length(row)
        if row(i) ~= col(i)
            if abs(row(i) - col(i)) - length(LabNet)/2 == 0
                % points are diametric, so draw a straight line
                u = [cos(theta(row(i)));sin(theta(row(i)))];
                v = [cos(theta(col(i)));sin(theta(col(i)))];
                %                 hline = line(...
                %                     [u(1);v(1)],...
                %                     [u(2);v(2)],...
                %                     'LineWidth', LG(i),...
                %                     'Color', (LG(i)/MaxLG)*edgeColors(i,:),...
                %                     'PickableParts','none')
                
                X = linspace(u(1),v(1),100);
                Y = linspace(u(2),v(2),100);
                Z = Y*0;
                [x,y,z]=tubeplot([X(:) Y(:) Z(:)]',LG(i)/350,20);
                XY = [x(1,:)' y(1,:)' z(1,:)'];
                t = cumsum(sqrt(sum((XY(1:end-1,1) -XY(2:end,1)).^2,2)));
                t =[t;max(t)];
                
                if directed
                    if invbool
                        if (row(i)<=idxs(1)) & (col(i)<=idxs(1))
                            colo = Val2colors(flipdim(t,1),cl_left);
                        elseif (row(i)>= idxs(1)+1) & (col(i)>= idxs(1)+1)
                            colo = Val2colors(flipdim(t,1),cl_right);
                        else
                            colo = Val2colors(flipdim(t,1),cl_inter);
                        end
                    else
                        
                        if (row(i) <= idxs(1)) & (col(i) <= idxs(1))
                            colo = Val2colors(t,cl_left);
                        elseif (row(i) >= idxs(1)+1) & (col(i) >= (idxs(1)+1))
                            colo = Val2colors(t,cl_right);
                        else
                            colo = Val2colors(t,cl_inter);
                        end

                        %colo = colorGradient([0 0 0],stColors(col(i),:),length(t));
                    end
                    colMat(:,:,1) = repmat(colo(:,1)' , [size(x,1) 1]);
                    colMat(:,:,2) = repmat(colo(:,2)' , [size(x,1) 1]);
                    colMat(:,:,3) = repmat(colo(:,3)' , [size(x,1) 1]);
                else
                    
                    colMat(:,:,1) = ones(size(x))*edgeColors(LabNet(row(i),col(i)),1);
                    colMat(:,:,2) = ones(size(x))*edgeColors(LabNet(row(i),col(i)),2);
                    colMat(:,:,3) = ones(size(x))*edgeColors(LabNet(row(i),col(i)),3);
                end
                
                
                htube = surf(x,y,z,colMat);
                clear colMat;
                %                 if bid
                %                     htube.FaceColor = Colors(i,:);
                %                 end
                htube.EdgeColor = 'none';
                htube.FaceAlpha = LG(i)/MaxLG;
                
                
            else % points are not diametric, so draw an arc
                u  = [cos(theta(row(i)));sin(theta(row(i)))];
                v  = [cos(theta(col(i)));sin(theta(col(i)))];
                x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
                y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
                r  = sqrt(x0^2 + y0^2 - 1);
                thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
                thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
                
                if u(1) >= 0 && v(1) >= 0
                    % ensure the arc is within the unit disk
                    thetat = [linspace(max(thetaLim),pi,50),...
                        linspace(-pi,min(thetaLim),50)].';
                else
                    thetat = linspace(thetaLim(1),thetaLim(2)).';
                end
                
                %                 hline = line(...
                %                     r*cos(thetat)+x0,...
                %                     r*sin(thetat)+y0,...
                %                     'LineWidth', LG(i),...
                %                     'Color', edgeColors(i,:),...
                %                     'PickableParts','none');
                %
                
                
                [x,y,z]=tubeplot([r*cos(thetat)+x0 r*sin(thetat)+y0 zeros(length(thetat),1)]',LG(i)/350,20);
                XY = [x(1,:)' y(1,:)' z(1,:)'];
                t = cumsum(sqrt(sum((XY(1:end-1,1) -XY(2:end,1)).^2,2)));
                t =[t;max(t)];
                
                if directed
                    if invbool
                        if (row(i) <= idxs(1)) & (col(i) <= idxs(1))
                            colo = Val2colors(flipdim(t,1),cl_left);
                        elseif (row(i) >= idxs(1)+1) & (col(i) >= idxs(1)+1)
                            colo = Val2colors(flipdim(t,1),cl_right);
                        else
                            colo = Val2colors(flipdim(t,1),cl_inter);
                        end
                    else
                        
                        if (row(i) <= idxs(1)) & (col(i) <= idxs(1))
                            colo = Val2colors(t,cl_left);
                        elseif (row(i) >= idxs(1)+1) & (col(i)>= idxs(1)+1)
                            colo = Val2colors(t,cl_right);
                        else
                            colo = Val2colors(t,cl_inter);
                        end
                        
                        %colo = colorGradient([0 0 0],stColors(col(i),:),length(t));
                    end
                    colMat(:,:,1) = repmat(colo(:,1)' , [size(x,1) 1]);
                    colMat(:,:,2) = repmat(colo(:,2)' , [size(x,1) 1]);
                    colMat(:,:,3) = repmat(colo(:,3)' , [size(x,1) 1]);
                    
                else
                    
                    colMat(:,:,1) = ones(size(x))*edgeColors(LabNet(row(i),col(i)),1);
                    colMat(:,:,2) = ones(size(x))*edgeColors(LabNet(row(i),col(i)),2);
                    colMat(:,:,3) = ones(size(x))*edgeColors(LabNet(row(i),col(i)),3);
                end
                
                htube = surf(x,y,z,colMat);
                clear colMat;
                %                 if bid
                %                     htube.FaceColor = Colors(i,:);
                %                 end
                htube.EdgeColor = 'none';
                htube.FaceAlpha = LG(i)/MaxLG;
            end
        end
    end
end

rmaxLim = param.rho(1);
Nstruct = size(CMatrix,1);
%  Constructing Bullseye

%% =========================== Drawing Data rings ======================= %
ColorTotal = [0 0 0];
rorig = rmaxLim;
if exist('CharFile','var')
    if ischar(CharFile)
        if exist(CharFile,'file')
            
            % The matrices are defined from left to right...graph is CCW
            data = textread(CharFile,'%f','delimiter',';');
            fid = fopen(CharFile,'rt');
            line1 = fgetl(fid);
            ind = strfind(line1,';');
            Nc = length(ind) + 1;
            Nr = size(data,1)/Nc;
            data = reshape(data,[Nc Nr])';
            %data = data./repmat(max(data),[size(data,1) 1]);
            
            data = data(1:size(CMatrix,1),:);
        end
    else
        data = CharFile;
    end
    Nrings = size(data,2);
    for i = 1: Nrings        
        
        rminLim = (rorig + i*param.intringwidth + (i-1)*param.ringwidth);
        rmaxLim = (rorig + i*param.intringwidth + i*param.ringwidth);
        param.rho = [rminLim rmaxLim];
        minplaneimg = min(min(data(:,i))); % find minimum first.
        scaledimg = (floor(((data(:,i) - minplaneimg) ./ ...
            (max(max(data(:,i))) - minplaneimg)) * 255)); % perform scaling
        
        % convert the image to a true color image with the jet colormap.
        ColorTotal = [ColorTotal;eval(colorstrings(i,:))];
        colorimg = squeeze(ind2rgb(scaledimg,eval(colorstrings(i,:))));
        
        param.tht   =[(theta(1)-ang/2)*180/pi (theta(end)+ang/2)*180/pi];
        %         [X,Y]=gen_xy(data(:,i),param);
        %         CDat = zeros(param.N*2,size(CMatrix,1),3);
        %         CDat(:,:,1)=repmat(colorimg(:,1)',[param.N*2 1]);
        %         CDat(:,:,2)=repmat(colorimg(:,2)',[param.N*2 1]);
        %         CDat(:,:,3)=repmat(colorimg(:,3)',[param.N*2 1]);
        %         h = patch(X,Y,repmat(data(:,i)',[param.N*2 1])*0,CDat);
        
        %           A = interp1(1:length(data)',data,1:length(data)/(length(data)*90):length(data)')';
        Npoints = 100;
        switch hemiVar
            
            case 1
                
                % Left Side
                %Ninterval = [1:Nstruct/2;];
                Ninterval = [1:idxs(1)];
                tempVar = data(Ninterval,i);
                tempColor = colorimg(Ninterval,:);
                Norig = length(tempVar);
                for k = 1:Norig-1
                    B = linspace(tempVar(k),tempVar(k+1),Npoints);
                    [col] = colorGradient(tempColor(k,:),tempColor(k+1,:),Npoints);
                    if k == 1
                        A = B(:);
                        Colors = col;
                    else
                        A = [A;B(2:end)'];
                        Colors = [Colors;col(2:end,:)];
                    end
                end
                param.tht   =[(theta(Ninterval(1)))*180/pi (theta(Ninterval(end)))*180/pi];
                [X,Y]=gen_xy(A,param);
                CDat = zeros(param.N*2,length(A),3);
                CDat(:,:,1)=repmat(Colors(:,1)',[param.N*2 1]);
                CDat(:,:,2)=repmat(Colors(:,2)',[param.N*2 1]);
                CDat(:,:,3)=repmat(Colors(:,3)',[param.N*2 1]);
                h = patch(X,Y,repmat(lin(A)',[param.N*2 1])*0,CDat, 'LineStyle', 'none');
                h.EdgeColor = 'none';
                
                
                % Right Side
                %Ninterval = [Nstruct/2+1:Nstruct];
                Ninterval = [idxs(1)+1:idxs(2)+idxs(1)];
                tempVar = data(Ninterval,i);
                tempColor = colorimg(Ninterval,:);
                Norig = length(tempVar);
                for k = 1:Norig-1
                    B = linspace(tempVar(k),tempVar(k+1),Npoints);
                    [col] = colorGradient(tempColor(k,:),tempColor(k+1,:),Npoints);
                    if k == 1
                        A = B(:);
                        Colors = col;
                    else
                        A = [A;B(2:end)'];
                        Colors = [Colors;col(2:end,:)];
                    end
                end
                param.tht   =[(theta(Ninterval(1)))*180/pi (theta(Ninterval(end)))*180/pi];
                [X,Y]=gen_xy(A,param);
                CDat = zeros(param.N*2,length(A),3);
                CDat(:,:,1)=repmat(Colors(:,1)',[param.N*2 1]);
                CDat(:,:,2)=repmat(Colors(:,2)',[param.N*2 1]);
                CDat(:,:,3)=repmat(Colors(:,3)',[param.N*2 1]);
                h = patch(X,Y,repmat(lin(A)',[param.N*2 1])*0,CDat, 'LineStyle', 'none');
                h.EdgeColor = 'none';
                
                
%                 % Frontal Right
%                 Nfront = 6;
%                 Ninterval = [1:Nfront] + Nstruct/2;
%                 
%                 
%                 [X,Y]=gen_xy(A,param);
%                 CDat = zeros(param.N*2,length(A),3);
%                 CDat(:,:,1)=repmat(colorimg(:,1)',[param.N*2 1]);
%                 CDat(:,:,2)=repmat(colorimg(:,2)',[param.N*2 1]);
%                 CDat(:,:,3)=repmat(colorimg(:,3)',[param.N*2 1]);
%                 h = patch(X,Y,repmat(lin(A)',[param.N*2 1])*0,CDat);
%                 h.EdgeColor = 'none';
%                 
%                 
%                 
%                 
%                 [X,Y,CDat] = Create_Ring(Ninterval, colorimg, param);
%                 h=patch(X,Y,repmat(lin(data(Ninterval,i)')',[param.N*2 1])*0,CDat);
%                 
                
            otherwise
                
                
                
                tempVar = [data(:,i); data(1,i)];
                tempColor = [colorimg; colorimg(1,:)];
                Norig = length(tempVar);
                for k = 1:Norig-1
                    B = linspace(tempVar(k),tempVar(k+1),Npoints);
                    [col] = colorGradient(tempColor(k,:),tempColor(k+1,:),Npoints);
                    if k == 1
                        A = B(:);
                        Colors = col;
                    else
                        A = [A;B(2:end)'];
                        Colors = [Colors;col(2:end,:)];
                    end
                end
                A = [A(end-Npoints/2:end);A(1:end-1-Npoints/2)];
                colorimg = [Colors(end-floor(Npoints/2):end,:);Colors(1:end-1-floor(Npoints/2),:)];
                [X,Y]=gen_xy(A,param);
                %           n = length(colorimg);
                %           X0 = linspace (1, n, length(A));
                %           s = interp1(data,colorimg,A);
                %           s1 = interp1(data,flipdim(colorimg,1),X0);
                
                %           colorimg = (s + s)/2;
                CDat = zeros(param.N*2,length(A),3);
                CDat(:,:,1)=repmat(colorimg(:,1)',[param.N*2 1]);
                CDat(:,:,2)=repmat(colorimg(:,2)',[param.N*2 1]);
                CDat(:,:,3)=repmat(colorimg(:,3)',[param.N*2 1]);
                h = patch(X,Y,repmat(lin(A)',[param.N*2 1])*0,CDat, 'LineStyle', 'none');
                h.EdgeColor = 'none';
        end
    end
end
%% ==================== End of Drawing Data rings ======================= %

%% ========================== Creating Concatenated Colormaps =========== %
if size(ColorTotal,1)>1
    ColorTotal(1,:) = [];
    colormap(ColorTotal);
    h = colorbar;
    Nticks = 2*size(data,2) + 1; % Number of Ticks
    StartTicks = linspace(0,1,Nticks);
    Temp = [1:2:length(StartTicks)];Temp(1) = [];Temp(end) = [];
    Temp2 = [StartTicks(Temp) + param.textspace StartTicks(Temp) - param.textspace];StartTicks(Temp) = [];
    Ticks = sort([StartTicks Temp2]);
    set(h, 'YTick', Ticks);
    lab = [min(data);mean(data);max(data)];
    labels = num2str(lab(:));
    set(h, 'YTickLabel', labels);
    set(get(h,'title'),'string','Multiple Colorbars');
    
end
%% ==================== End of Creating Colormap ======================== %

%% ================= Drawing Rings for Structures ======================= %
param.rho = [rmaxLim + param.intringwidth rmaxLim + param.intringwidth + 0.1];
rmaxLim = rmaxLim+param.intringwidth+0.1;
Ncolor = size(stColors,1);
Nstruct = size(CMatrix,1);
re = floor(size(CMatrix,1)/Ncolor);
stColors = repmat(stColors,[re+1 1]);
stColors = stColors(1:Nstruct,:);

% Unknown parcellation

switch hemiVar
    
    case 1
        
        % Frontal Left
        %Nleft = Nstruct/2;
        Nleft = idxs(1);
        Ninterval = [1:Nleft];
        param.tht   =[(theta(Ninterval(1))-ang/2)*180/pi (theta(Ninterval(end))+ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, stColors, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat, 'LineStyle', 'none');
        
        
        %Ninterval = [1:Nleft] + Nstruct/2;
        Ninterval = [idxs(1)+1:idxs(2)+idxs(1)];
        param.tht   =[(theta(Ninterval(1))+ang/2)*180/pi (theta(Ninterval(end))-ang/2)*180/pi];
        [X,Y,CDat] = Create_Ring(Ninterval, stColors, param);
        h=patch(X,Y,CDat(:,:,1)*0,CDat, 'LineStyle', 'none');
        %theta2 = linspace(pi/2-ang,-pi/2+ang,(idxs(2)))'; %theta2(end) = [];   
        

        

    otherwise
        param.tht   =[(theta(1)-ang/2)*180/pi (theta(end)+ang/2)*180/pi];
        [X,Y]=gen_xy([1:Nstruct]',param);
        CDat = zeros(param.N*2,Nstruct,3);
        CDat(:,:,1)=repmat(stColors(:,1)',[param.N*2 1]);
        CDat(:,:,2)=repmat(stColors(:,2)',[param.N*2 1]);
        CDat(:,:,3)=repmat(stColors(:,3)',[param.N*2 1]);
        h = patch(X,Y,squeeze(CDat(:,:,1))*0,CDat,'LineStyle', 'none');
end

%% ============== End of Drawing Rings for Structures ====================%
xy_text = rmaxLim.* [cos(theta) sin(theta)];
% 
% switch hemiVar
%     case 1
        ipos       = xy_text(:,1) > 0;
        s.t        = zeros(Nstruct,1);
        s.t( ipos) = text(xy_text(ipos,1).*1.05, xy_text(ipos,2).*1.05, stNames(ipos), 'FontSize',16,'FontName','Helvetica', 'interpreter', 'latex');
        set(s.t( ipos), {'Rotation'},num2cell(theta(ipos)*180/pi));
        s.t( ~ipos) = text(xy_text(~ipos,1).*1.05, xy_text(~ipos,2).*1.05, stNames(~ipos), 'FontSize',16,'FontName','Helvetica', 'interpreter', 'latex');
        set(s.t( ~ipos), {'Rotation'},num2cell(theta(~ipos)*180/pi-180),'Horiz','right');
%     otherwise
%         h = text(xy_text(:,1).*1.05, xy_text(:,2).*1.05, stNames, 'FontSize',12);
%         set(h, {'Rotation'},num2cell(theta*180/pi));
% end

camlight
return;

function [X,Y,CDat] = Create_Ring(Ninterval, Colors, param);
%
% Syntax :
% [X,Y,CDat] = Create_Ring(Ninterval,theta, Colors, param);
%
% This function creates the coordinates and color data for bullseye
% plotting.
%
% Input Parameters:
%   Ninterval      : Structures interval
%   theta          : Theta angles (radians)
%   Colors         : Colors
%   param          : Parameters
%
%
% Output Parameters:
%   X              : X coordinates
%   Y              : Y coordinates
%   CDat           : Color Data
%
% Related references:
%
%
% See also: Plot_Circular_Graph
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 12th 2013
% Version $1.0

Nlength = length(Ninterval);
CDat = zeros(param.N*2,Nlength,3);
[X,Y]=gen_xy(Ninterval',param);
CDat(:,:,1)=repmat(Colors(Ninterval,1)',[param.N*2 1]);
CDat(:,:,2)=repmat(Colors(Ninterval,2)',[param.N*2 1]);
CDat(:,:,3)=repmat(Colors(Ninterval,3)',[param.N*2 1]);
return


%% Bullseye internal functions

function [X,Y]=gen_xy(data,param);
ind=1;
for j=1:size(data,1)
    for k=1:size(data,2)
        dtht=diff(param.tht)/size(data,1);
        tht1=dtht*(j-1)+param.tht(1)-param.tht0;
        tht2=dtht* j   +param.tht(1)-param.tht0;
        
        dr=diff(param.rho)/size(data,2);
        rho1=dr*(k-1)+param.rho(1);
        rho2=dr* k   +param.rho(1);
        
        ang = linspace(tht1/180*pi,tht2/180*pi,param.N);
        [arc1.x,arc1.y]=pol2cart(ang,rho1);
        [arc2.x,arc2.y]=pol2cart(ang,rho2);
        X(:,ind)=[arc1.x arc2.x(end:-1:1)];
        Y(:,ind)=[arc1.y arc2.y(end:-1:1)];
        ind=ind+1;
    end
end
return

function t=bull_labels(param);

% Labels are added CCW...this alone determines that the first ROW is INFEROLATERAL etc...
% label{1}='inferoseptal';
% label{2}='inferior';
% label{3}='inferolateral';
% label{4}='anterolateral';
% label{5}='anterior';
% label{6}='anteroseptal';
label{1}='anterior';
label{2}='anteroseptal';
label{3}='inferoseptal';
label{4}='inferior';
label{5}='inferolateral';
label{6}='anterolateral';

for k=1:length(label)
    %   ang=(60*(k-1)-param.tht0-30)*(pi/180);
    ang=(60*(k-1)-param.tht0-60)*(pi/180);
    %   [x,y]=pol2cart(ang,1.08*max(param.rho));
    [x,y]=pol2cart(ang,0.875*max(param.rho));
    t(k)=text(x,y,label{k});
    if rem((ang*180/pi),360)>=180
        set(t(k),'rotation',ang*(180/pi)+90,'HorizontalAlignment','Center','FontSize',14,'FontName','Courier','FontWeight','Bold');
    elseif (ang*180/pi)==0 | (ang*180/pi)==360
        set(t(k),'rotation',ang*(180/pi)+90,'HorizontalAlignment','Center','FontSize',14,'FontName','Courier','FontWeight','Bold');
    else
        set(t(k),'rotation',ang*(180/pi)-90,'HorizontalAlignment','Center','FontSize',14,'FontName','Courier','FontWeight','Bold');
    end
end

return

%  This function linearizes a matrix (m) of any dimension (eg M=m(:)).
%  If an index vector (ind) is given then the the ind entries of m(:) are
%  returned.
%
% SYNTAX: m=lin(m);
%         m=lin(m,ind);
%
% DBE 12/22/03

function m=lin(m,ind);

m=m(:);

if nargin==2
    m=m(ind);
end

function [GMcodes,Names, Colors] = Brain_GM_codes(atlastype);
%
% Syntax :
% [GMcodes,Names, Colors] = Brain_codes(atlastype);
%
% This function extract all gray matter codes from a specified atlas.
%
% Input Parameters:
%   atlastype     : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
%   Names         : Structure Names
%   Colors        : Structure Colors
%
% Related references:
%
%
% See also:
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0
GMcodes = '';
Names = '';
cont = 0;
Colors = [0 0 0];
switch atlastype
    case 'aparc+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [7 8 10:13 17:18  26 46 47 49:54 58 1001:1003 1005:1035 2001:2003 2005:2035];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2009s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [7 8 10:13 17:18  26 46 47 49:54 58 11101:11175 12101:12175];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2005s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [7 8 10:13 17:18  26 46 47 49:54 58 1102:1181 2102:2181];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'yeoh7networks'
        
    case 'yeoh17networks'
        
    case 'parckmeans'
        
    case 'ibaspm116'
        txtfile = which('atlas116.cod');
        col = 255*[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        re = floor(length(Codes)/Ncolor); col = repmat(col,[re+1 1]);
        %         Codes = str2num(char(a));
        Names = char(b);
        Colors = col(1:length(Codes),:);
        GMcodes = Codes;
    case 'ibaspm71'
        txtfile = which('atlas71.cod');
        col = 255*[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        re = floor(length(Codes)/Ncolor); col = repmat(col,[re+1 1]);
        Colors = col(1:length(Codes),:);
        GMcodes = Codes;
    otherwise
        
end
return

function NETWORKS = Label_Networks(conn, threshold);
%
% Syntax :
%  Networks = Label_Networks(conn, threshold);
%
% This script labels the unconnected networks in a connectivity matrix
% provided they are higher than a certain threshold
%
% Input Parameters:
%   conn              :  Connectivity Matrix
%   threshold         :  Threshold of the minimum intensity in the
%                        connectivity matrix
%
% Output Parameters:
%   NETWORKS          :  Matrix that contains a label in every element
%                        in conn matrix according to the network they belong
%
% Related references:
%
%
% See also: graphshortestpath
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% Jul 27th 2012
% Version $1.0

% ======================= Main Program ================================== %
% conn = triu(conn);
%  conn_t = conn > threshold;
%  dims = size(conn);
%  NETWORKS = zeros(dims);
%  indices = find(conn_t>0);
%  num = 1;
%
%  [ix,iy] = ind2sub(size(conn_t),indices);
%  index = [ix iy];
%  Id = ix*0;
%  while sum(logical(Id)) ~= length(Id)
%      idf = index(1,:);
%      while ~isempty(idf)
%          indcon = find(sum(ismember(index,idf)') > 0);
%          idf = index(indcon,:);
%          Id(indcon) = num;
%          index(indcon,:) = [];
%          idf = unique(idf(:));
%      end
%      index(find(Id == num),:) = [];
%      num = num +1;
%  end
%      indices = find(conn_t>0);
% %
% %     elems = [];
%
%     while(length(elems)>0)
%         [X, Y] = ind2sub(dims, elems(1));
%         elems(1)=[];
%         NETWORKS(X,Y)=num;
%         NETWORKS(Y,X)=num;
%         conn_t(X,Y) = 0;
%         conn_t(Y,X) = 0;
%         a = find(conn_t(Y,:)>0);
%         for i = 1:length(a)
%             if (a(i)~=X)
%                 elems = [elems sub2ind(dims,Y,a(i))];
%             end
%         end
%         a = find(conn_t(X,:)>0);
%         for i = 1:length(a)
%             if (a(i)~=Y)
%                 elems = [elems sub2ind(dims,a(i),X)];
%             end
%         end
%
%     end
%     num = num+1;
%
%  end




conn_t = conn > threshold;
dims = size(conn);
NETWORKS = zeros(dims);
indices = find(conn_t>0);
num = 1;


while(sum(conn_t(:))>0)
    
    indices = find(conn_t>0);
    
    elems = [indices(1)];
    
    while(length(elems)>0)
        [X, Y] = ind2sub(dims, elems(1));
        elems(1)=[];
        NETWORKS(X,Y)=num;
        NETWORKS(Y,X)=num;
        conn_t(X,Y) = 0;
        conn_t(Y,X) = 0;
        a = find(conn_t(Y,:)>0);
        for i = 1:length(a)
            if (a(i)~=X)
                elems = [elems sub2ind(dims,Y,a(i))];
            end
        end
        a = find(conn_t(X,:)>0);
        for i = 1:length(a)
            if (a(i)~=Y)
                elems = [elems sub2ind(dims,a(i),X)];
            end
        end
        
    end
    num = num+1;
    
end

function [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Syntax :
% [labid] = Recur_Corr(Surf,lab,labid,cont);
%
% Recursive function for labelling correction.
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Plot_oversurf
% Exp_Surf
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% December 1st 2007
% Version $1.0

ind = find(Surf.Is ==lab&(labid ==0));
A = unique(Surf.SurfData.faces(Surf.Tri(ind(1),3:2+Surf.Tri(ind(1),2)),:));A = A(:);
indt = find(Surf.Is(A)~= lab);
A(indt) = [];A(A == ind(1))= [];
T = unique([ind(1); A]);
labid(T) = cont;
An = rand(size(A));
while sum(A)~=sum(An)
    An = A;
    Neib = Surf.Tri(A,3:end); Neib = unique(nonzeros(Neib(:)));
    A = unique(Surf.SurfData.faces(Neib,:));A = A(:);
    indt = find(Surf.Is(A)~= lab);
    A(indt) = [];
    labid(A) = cont;
    T =unique([T;A]);
end
indn = find(Surf.Is ==lab&(labid ==0));
if ~isempty(indn)
    cont = cont+1;
    [labid] = Recur_Corr(Surf,lab,labid,cont);
else
    return;
end
if ~isempty(A)
    for i = 1:size(A,2)
        TA = unique(Surf.SurfData.faces(Surf.Tri(A(i),3:2+Surf.Tri(A(i),2)),:));
        indt = find(Surf.Is(TA(i))~= lab);
        TA(indt) = [] ;
        T = [T; A];
    end
else
    return;
end
return

function AcronSt  = Names_2_Acron(atype,Names);
%
% Syntax :
% AcronSt  = Names_2_Acron(atype,Names);
%
% This function creates anocryms for structures name according with a
% specified atlas type.
%
% Input Parameters:
%   atype          : Atlas type
%   Names          : Structure Names
%
% Output Parameters:
%     AcronSt      : Acronym names
%
% Related references:
%
%
% See also: Plot_Circular_Graph FreeSurf_2_Circos
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 12th 2013
% Version $1.0

switch atype
    case 'aparc+aseg'
        AcronSt = '';
        for i = 1:size(Names,1);
            
            % Left Hemisphere
            if strcmp(deblank(Names(i,:)),'Left-Thalamus')%|strcmp(deblank(Names(i,:)),'Right-Thalamus')
                AcronSt = strvcat(AcronSt,'LH-Thalam');
            elseif strcmp(deblank(Names(i,:)),'Left-Thalamus-Proper')%|strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')
                AcronSt = strvcat(AcronSt,'LH-ThalProp');
            elseif strcmp(deblank(Names(i,:)),'Left-Caudate')%|strcmp(deblank(Names(i,:)),'Right-Caudate')
                AcronSt = strvcat(AcronSt,'LH-Caudt');
            elseif strcmp(deblank(Names(i,:)),'Left-Putamen')%|strcmp(deblank(Names(i,:)),'Right-Putamen')
                AcronSt = strvcat(AcronSt,'LH-Putam');
            elseif strcmp(deblank(Names(i,:)),'Left-Pallidum')%|strcmp(deblank(Names(i,:)),'Right-Pallidum')
                AcronSt = strvcat(AcronSt,'LH-Pallid');
            elseif strcmp(deblank(Names(i,:)),'Left-Hippocampus')%|strcmp(deblank(Names(i,:)),'Right-Hippocampus')
                AcronSt = strvcat(AcronSt,'LH-Hippoc');
            elseif strcmp(deblank(Names(i,:)),'Left-Amygdala')%|strcmp(deblank(Names(i,:)),'Right-Amygdala')
                AcronSt = strvcat(AcronSt,'LH-Amygd');
            elseif strcmp(deblank(Names(i,:)),'Left-Accumbens-area')%|strcmp(deblank(Names(i,:)),'Right-Accumbens-area')
                AcronSt = strvcat(AcronSt,'LH-AccumbArea');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-bankssts')%|strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')
                AcronSt = strvcat(AcronSt,'LH-BanksStS');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-caudalanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')
                AcronSt = strvcat(AcronSt,'LH-CaudAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-caudalmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')
                AcronSt = strvcat(AcronSt,'LH-CaudMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-cuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')
                AcronSt = strvcat(AcronSt,'LH-Cuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-entorhinal')%|strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')
                AcronSt = strvcat(AcronSt,'LH-Entorh');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-fusiform')%|strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')
                AcronSt = strvcat(AcronSt,'LH-Fusif');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-inferiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')
                AcronSt = strvcat(AcronSt,'LH-InfPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-inferiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')
                AcronSt = strvcat(AcronSt,'LH-InfTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-isthmuscingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')
                AcronSt = strvcat(AcronSt,'LH-IsthCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lateraloccipital')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')
                AcronSt = strvcat(AcronSt,'LH-LatOccip');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lateralorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')
                AcronSt = strvcat(AcronSt,'LH-LatOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-lingual')%|strcmp(deblank(Names(i,:)),'ctx-rh-lingual')
                AcronSt = strvcat(AcronSt,'LH-Lingual');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-medialorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')
                AcronSt = strvcat(AcronSt,'LH-MedOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-middletemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')
                AcronSt = strvcat(AcronSt,'LH-MidTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parahippocampal')%|strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')
                AcronSt = strvcat(AcronSt,'LH-Parahip');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-paracentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')
                AcronSt = strvcat(AcronSt,'LH-Paracent');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parsopercularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')
                AcronSt = strvcat(AcronSt,'LH-ParsoPerc');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parsorbitalis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')
                AcronSt = strvcat(AcronSt,'LH-ParsOrb');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-parstriangularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')
                AcronSt = strvcat(AcronSt,'LH-ParsTriang');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-pericalcarine')%|strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')
                AcronSt = strvcat(AcronSt,'LH-PeriCalc');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-postcentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')
                AcronSt = strvcat(AcronSt,'LH-PostCent');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-posteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')
                AcronSt = strvcat(AcronSt,'LH-PostCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-precentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-precentral')
                AcronSt = strvcat(AcronSt,'LH-Precent');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-precuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')
                AcronSt = strvcat(AcronSt,'LH-Precuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-rostralanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')
                AcronSt = strvcat(AcronSt,'LH-RostAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-rostralmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')
                AcronSt = strvcat(AcronSt,'LH-RostMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiorfrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')
                AcronSt = strvcat(AcronSt,'LH-SupFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')
                AcronSt = strvcat(AcronSt,'LH-SupPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-superiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')
                AcronSt = strvcat(AcronSt,'LH-SupTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-supramarginal')%|strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')
                AcronSt = strvcat(AcronSt,'LH-Supramarg');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-frontalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')
                AcronSt = strvcat(AcronSt,'LH-FrontPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-temporalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')
                AcronSt = strvcat(AcronSt,'LH-TempPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-transversetemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')
                AcronSt = strvcat(AcronSt,'LH-TransvTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-lh-insula')%|strcmp(deblank(Names(i,:)),'ctx-rh-insula')
                AcronSt = strvcat(AcronSt,'LH-Insula');
                
                % Right Hemisphere
                
            elseif strcmp(deblank(Names(i,:)),'Right-Thalamus')%|strcmp(deblank(Names(i,:)),'Right-Thalamus')
                AcronSt = strvcat(AcronSt,'RH-Thalam');
            elseif strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')%|strcmp(deblank(Names(i,:)),'Right-Thalamus-Proper')
                AcronSt = strvcat(AcronSt,'RH-ThalProp');
            elseif strcmp(deblank(Names(i,:)),'Right-Caudate')%|strcmp(deblank(Names(i,:)),'Right-Caudate')
                AcronSt = strvcat(AcronSt,'RH-Caudt');
            elseif strcmp(deblank(Names(i,:)),'Right-Putamen')%|strcmp(deblank(Names(i,:)),'Right-Putamen')
                AcronSt = strvcat(AcronSt,'RH-Putam');
            elseif strcmp(deblank(Names(i,:)),'Right-Pallidum')%|strcmp(deblank(Names(i,:)),'Right-Pallidum')
                AcronSt = strvcat(AcronSt,'RH-Pallid');
            elseif strcmp(deblank(Names(i,:)),'Right-Hippocampus')%|strcmp(deblank(Names(i,:)),'Right-Hippocampus')
                AcronSt = strvcat(AcronSt,'RH-Hippoc');
            elseif strcmp(deblank(Names(i,:)),'Right-Amygdala')%|strcmp(deblank(Names(i,:)),'Right-Amygdala')
                AcronSt = strvcat(AcronSt,'RH-Amygd');
            elseif strcmp(deblank(Names(i,:)),'Right-Accumbens-area')%|strcmp(deblank(Names(i,:)),'Right-Accumbens-area')
                AcronSt = strvcat(AcronSt,'RH-AccumbArea');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')%|strcmp(deblank(Names(i,:)),'ctx-rh-bankssts')
                AcronSt = strvcat(AcronSt,'RH-BanksStS');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalanteriorcingulate')
                AcronSt = strvcat(AcronSt,'RH-CaudAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-caudalmiddlefrontal')
                AcronSt = strvcat(AcronSt,'RH-CaudMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-cuneus')
                AcronSt = strvcat(AcronSt,'RH-Cuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')%|strcmp(deblank(Names(i,:)),'ctx-rh-entorhinal')
                AcronSt = strvcat(AcronSt,'RH-Entorh');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')%|strcmp(deblank(Names(i,:)),'ctx-rh-fusiform')
                AcronSt = strvcat(AcronSt,'RH-Fusif');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiorparietal')
                AcronSt = strvcat(AcronSt,'RH-InfPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-inferiortemporal')
                AcronSt = strvcat(AcronSt,'RH-InfTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-isthmuscingulate')
                AcronSt = strvcat(AcronSt,'RH-IsthCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateraloccipital')
                AcronSt = strvcat(AcronSt,'RH-LatOccip');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-lateralorbitofrontal')
                AcronSt = strvcat(AcronSt,'RH-LatOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-lingual')%|strcmp(deblank(Names(i,:)),'ctx-rh-lingual')
                AcronSt = strvcat(AcronSt,'RH-Lingual');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-medialorbitofrontal')
                AcronSt = strvcat(AcronSt,'RH-MedOrbFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-middletemporal')
                AcronSt = strvcat(AcronSt,'RH-MidTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')%|strcmp(deblank(Names(i,:)),'ctx-rh-parahippocampal')
                AcronSt = strvcat(AcronSt,'RH-Parahip');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-paracentral')
                AcronSt = strvcat(AcronSt,'RH-Paracent');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsopercularis')
                AcronSt = strvcat(AcronSt,'RH-ParsoPerc');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parsorbitalis')
                AcronSt = strvcat(AcronSt,'RH-ParsOrb');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')%|strcmp(deblank(Names(i,:)),'ctx-rh-parstriangularis')
                AcronSt = strvcat(AcronSt,'RH-ParsTriang');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')%|strcmp(deblank(Names(i,:)),'ctx-rh-pericalcarine')
                AcronSt = strvcat(AcronSt,'RH-PeriCalc');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-postcentral')
                AcronSt = strvcat(AcronSt,'RH-PostCent');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-posteriorcingulate')
                AcronSt = strvcat(AcronSt,'RH-PostCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-precentral')%|strcmp(deblank(Names(i,:)),'ctx-rh-precentral')
                AcronSt = strvcat(AcronSt,'RH-Precent');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')%|strcmp(deblank(Names(i,:)),'ctx-rh-precuneus')
                AcronSt = strvcat(AcronSt,'RH-Precuneus');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralanteriorcingulate')
                AcronSt = strvcat(AcronSt,'RH-RostAntCing');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-rostralmiddlefrontal')
                AcronSt = strvcat(AcronSt,'RH-RostMidFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorfrontal')
                AcronSt = strvcat(AcronSt,'RH-SupFront');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiorparietal')
                AcronSt = strvcat(AcronSt,'RH-SupPariet');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-superiortemporal')
                AcronSt = strvcat(AcronSt,'RH-SupTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')%|strcmp(deblank(Names(i,:)),'ctx-rh-supramarginal')
                AcronSt = strvcat(AcronSt,'RH-Supramarg');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-frontalpole')
                AcronSt = strvcat(AcronSt,'RH-FrontPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')%|strcmp(deblank(Names(i,:)),'ctx-rh-temporalpole')
                AcronSt = strvcat(AcronSt,'RH-TempPole');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')%|strcmp(deblank(Names(i,:)),'ctx-rh-transversetemporal')
                AcronSt = strvcat(AcronSt,'RH-TransvTemp');
            elseif strcmp(deblank(Names(i,:)),'ctx-rh-insula')%|strcmp(deblank(Names(i,:)),'ctx-rh-insula')
                AcronSt = strvcat(AcronSt,'RH-Insula');
            end
        end
end
return;

function [GMcodes,Names, Colors] = Gray_Matter_codes(atlastype);
%
% Syntax :
% [GMcodes] = Gray_Matter_codes(atlastype);
%
% This function extract cortical gray matter codes from a specified atlas.
%
% Input Parameters:
%   atlastype     : Atlas type.
%
%
% Output Parameters:
%   GMcodes       : Gray matter codes.
%   Names         : Structure Names
%   Colors        : Structure Colors
%
% Related references:
%
%
% See also:
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0
GMcodes = '';
Names = '';
cont = 0;
Colors = [0 0 0];
col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
Ncolor = size(col,1);
switch atlastype
    case 'aparc+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 1001:1003 1005:1035 2001:2003 2005:2035];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2009s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 11101:11175 12101:12175];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'a2005s+aseg'
        [a,temp] = system('echo $FREESURFER_HOME');
        txtfile =  strcat(temp,filesep,'FreeSurferColorLUT.txt');
        fid = fopen(txtfile);
        cont = 0;lines = '';
        while 1
            cont = cont + 1;
            line = fgetl(fid);
            if ~ischar(line),   break,   end
            lines = strvcat(lines,line);
        end
        ind = find(sum(isspace(lines'))==size(lines,2));lines(ind,:) = [];
        Codes = '';
        Names = '';
        cont = 0;
        Colors = [0 0 0];
        for i = 1:size(lines,1)
            if ~strcmp(lower(deblank(lines(i,1))),'#');
                cont = cont + 1;
                temp = strread(deblank(lines(i,:)),'%s','delimiter',' ');
                code = temp{1};
                name = temp{2};
                Codes = strvcat(Codes,code);
                Names = strvcat(Names,name);
                Colors(cont,:) = [str2num(temp{3}) str2num(temp{4}) str2num(temp{5})];
            end
        end
        fclose(fid);
        GMcodes = [10:13 17:18 26 49:54 58 1102:1181 2102:2181];
        index = ismember(str2num(Codes),GMcodes);
        Names = Names(index,:);
        Colors = Colors(index,:);
    case 'yeoh7networks'
        
    case 'yeoh17networks'
        
    case 'parckmeans'
        GMcodes = [10:13 17:18 26 49:54 58 11101:11175 12101:12175];
        re = floor(90/Ncolor); col = repmat(col,[re+1 1]);
    case 'ibaspm116'
        txtfile = which('atlas116.cod');
        re = floor(90/Ncolor); col = repmat(col,[re+1 1]);
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1:90];
        index = ismember(Codes,GMcodes);
        Names = Names(index,:);
        Colors = col(1:90,:);
    case 'ibaspm71'
        txtfile = which('atlas71.cod');
        col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
        Ncolor = size(col,1);
        [Codes,b] = textread(txtfile,'%u%s','delimiter','=');
        %         Codes = str2num(char(a));
        Names = char(b);
        GMcodes = [1     2     4     5     6     7     9    10    11    12    14    15    16    18    19    20 ...
            23    25    26    27    32    33    36    37    38    39    41    50    52    53    54    56    60 ...
            61    62    63    64    67    69    70    72    74    75    76    80    85    88    90    97    98 ...
            99   101   102   108   110   112   114   119   125   130   132   140   145   154   159   164   165 ...
            175   196   203   251];
        index = ismember(Codes,GMcodes);
        re = floor(length(index)/Ncolor); col = repmat(col,[re+1 1]);
        Names = Names(index,:);
        Colors = col(1:length(index),:);
end

return


function  Connect = Read_Connectivity_Matrix(ConnectFilename,vis,cl);
%
% Syntax :
%  Connect = Read_Connectivity_Matrix(ConnectFilename,vis,cl);
%
% This script reads the connectivity Matrix from a connectivity matrix filename.
%
% Input Parameters:
%   ConnectFilename   :  Text Connectivity Matrix Filename
%      vis            :  Boolean Variable to plot, or not, the Connectivity
%                        Matrix
%        cl           :  Colormap
%
% Output Parameters:
%   Connect           :  Connectivity Matrix
%
% Related references:
%
%
% See also:
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% May 18th 2012
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin ==0
    error('Please select a correct connectivity matrix');
    return;
end
if nargin <2
    vis = 1;
    cl = 'jet';
end
if nargin ==2
    cl = 'jet';
end
%=====================End of Checking Input Parameters====================%
%==========================  Main Program   ==============================%
if ~exist(ConnectFilename,'file')
    error('The specified Matrix file does not exist');
    return;
end
fid  = fopen(ConnectFilename,'rt');
Nstruct = str2num(char(textread(ConnectFilename,'%s',1,'headerlines',8)));
Snames = textread(ConnectFilename,'%s',Nstruct*2,'headerlines',9);
temp = char(Snames); names= temp(2:2:end,:);codes= temp(1:2:end,:);
temp = textread(ConnectFilename,'%f',Nstruct*Nstruct,'headerlines',Nstruct+11);
Connect.Matrix = reshape(temp,[Nstruct Nstruct]);
Connect.StructCodes = str2num(deblank(codes));
Connect.StructNames = deblank(names);
if vis
    [pth,nm,ext] = fileparts(ConnectFilename);
    colordef black;figure;set(gcf,'Color',[0 0 0]);
    imagesc(Connect.Matrix);hold on;plot([1:Nstruct],[1:Nstruct],'xw');colormap(cl);
    set(gca,'YTick',[1:Nstruct]);set(gca,'XTick',[1:Nstruct]); set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',Connect.StructNames);
    title([nm ' ===> ' num2str(Nstruct) ' Structures']);
    
    h = colorbar;
    if ~strcmp(cl,'spectral')
        colormap(cl);
    else
        colormap(spectral(64));
    end
    range = max(Connect.Matrix(:))-min(Connect.Matrix(:));
    values =  min(Connect.Matrix(:)):range/10:max(Connect.Matrix(:));
    set(h,'YTickLabel',num2str(values'));
    
end
%=============================End=========================================%
return;
