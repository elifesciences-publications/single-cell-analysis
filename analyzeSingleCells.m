%% Analyze single cell biofilm parameters as heatmaps
%
% FORMAT analyzeSingleCells(cells, varargin)
% 
% cells    - [required] structure containing the fields "data", "measurements", 
%               "scaling_dxyz" and "time" (in s)
%               
%               - data: 3D cell data in the format that the Matlab function 
%                 "bwconncomp" (Image Processing Toolbox) would return with
%                 fields "Connectivity", "ImageSize", "NumObjects" and
%                 "PixelIdxList" (see bwconncomp function documentation for more
%                 information). "labelmatrix(cells(i).data)" will create a
%                 3D matrix with individual cells labelled
%               - measurements: structure containing single cell 
%                 measurements. The line index is referring to the object 
%                 index of the "data" structure. 
%
%                   the example data ("example_biofilm_R234A.mat") contains the following fields: 
%                     - Centroid: Centroid of each cell as obtained by the
%                       Matlab function "regionprops" (Image Processing
%                       Toolbox) [voxels]
%                     - BoundingBox: BoundingBox of each cell as obtained by the
%                       Matlab function "regionprops" (Image Processing 
%                       Toolbox) [voxels]
%                     - Volume: Volume of each cell [µm^3]
%                     - distanceToCM_xyz: Distance to center of mass [µm]
%                     - distanceToCM_xy: Distance to center of mass projected 
%                       onto the surface of the subtrate [µm]
%                     - distanceToNearestNeighbor: Nearest neighbor
%                       distance (centroid to centroid) [µm]
%                     - width/length/height: Width/length/height of each cell
%                       obtained by principal component anaylsis (PCA) [µm]
%                     - evecs: Eigenvectores of each cell obtained by PCA
%                     - alpha_flow: Angle between the cell's major axis and the
%                       direction of the external growth medium medium flow
%                       field [rad]
%                     - alpha_z: Angle between the cell's major axis and z
%                       [rad]
%                     - alpha_radius: Angle between the cell's major axis
%                       and the vector pointing into the cells direction from
%                       the biofilm's center of mass projected onto the
%                       subtrate
%                     - nematicOrderParameter_3: Nematic order parameter in
%                       a local neighborhood of 3 µm
%                     - distanceToSurface: Distance to biofilm surface [µm]
%
%               - scaling_dxyz: voxel spacing in µm/voxel
%               - time_interval: duration between frame i and i-1 [s]
% 
% fieldY   - [optional] measurement field to plot on the y-axis (for available 
%            fields of the example data see above) (default: 'z')
% fieldZ   - [optional] measurement field to plot as colored tiles (for available
%            fields of the example data see above) (default:
%            'nematicOrderParameter_3')
% rangeY   - [optional] range of Y (default: [min max])
% rangeZ   - [optional] range of Z/colormap range (default: [min max])
% nBinsY   - [optional] number of bins along Y (default: 30)
% averagingFcn   - [optional] function to average data points for each heatmap tile 
%                  (default: 'mean', other options: 'median', 
%                  'sum', 'min', 'max')
%
% _________________________________________________________________________
%% Copyright (C) 2017 Max Planck Institute for Terrestrial Microbiology,
% Feel free to contact Raimo Hartmann (raimo.hartmann@gmail.com)
% for more information
% _________________________________________________________________________
%
%% Example (create heatmap displaying the nematic order vs. biofilm height vs. time:
%   load('example_biofilm_R234A.mat')
%   analyzeSingleCells(cells, 'fieldY', 'z', 'fieldZ', 'nematicOrderParameter_3', 'rangeY', [0 20], 'rangeZ', [-0.5 1], 'nBinsY', 30)
%

%%
function analyzeSingleCells(cells, varargin)
% Handle inputs
p = inputParser;

addRequired(p,'cells');
addParameter(p,'fieldY', 'z');
addParameter(p,'fieldZ', 'nematicOrderParameter_3');
addParameter(p,'rangeY', []);
addParameter(p,'rangeZ', []);
addParameter(p,'nBinsY', 30);
addParameter(p,'averagingFcn', 'mean');

parse(p,cells,varargin{:});

fieldY = p.Results.fieldY;
fieldZ = p.Results.fieldZ;
rangeY = p.Results.rangeY;
rangeZ = p.Results.rangeZ;
nBinsY = p.Results.nBinsY;
averagingFcn = p.Results.averagingFcn;

clear p;

% Check inputs
if isempty(find(strcmp(fieldnames(cells), 'data'), 1))
    error('Input structure is missing the field "data"');
end
if isempty(find(strcmp(fieldnames(cells), 'measurements'), 1))
    error('Input structure is missing the field "measurements"');
end
if isempty(find(strcmp(fieldnames(cells), 'scaling_dxyz'), 1))
    error('Input structure is missing the field "scaling_dxyz"');
end
if ~strcmp(fieldY, 'z') && ~find(strcmp(fieldY, fieldnames(cells(1).measurements))) || ~isempty(find(strcmp(fieldY, {'Centroid', 'BoundingBox', 'evecs'}), 1))
    error('FieldY is not containing a valid single cell measurement for plotting');
end
if ~strcmp(fieldZ, 'z') && ~find(strcmp(fieldZ, fieldnames(cells(1).measurements))) || ~isempty(find(strcmp(fieldY, {'Centroid', 'BoundingBox', 'evecs'}), 1)) 
    error('FieldZ is not containing a valid single cell measurement for plotting');
end

% Prepare data

% Remove z-offset
zOffset = zeros(1, numel(cells));
for i = 1:numel(cells)
    try
        centroids = [cells(i).measurements.Centroid];
        centroids_z = centroids(3:3:end);
        zOffset(i) = min(centroids_z)*cells(i).scaling_dxyz;
    catch
        zOffset(i) = NaN;
    end
end
zOffset = nanmedian(zOffset);

% x-axis
time_intervals = [cells.time_interval]/60/60;
time = cumsum(time_intervals); % [h]
rangeX = [time(1) time(end)+time_intervals(end)];
binsX = [time(find(time>=rangeX(1) & time<=rangeX(2))) time(find(time<=rangeX(2), 1, 'last'))+time_intervals(find(time<=rangeX(2), 1, 'last'))];

% y-axis
if isempty(rangeZ)
    minY = inf(numel(numel(cells), 1));
    maxY = -inf(numel(numel(cells), 1));
    for i = 1:numel(cells)
        switch fieldY
            case 'z'
                data_temp = [cells(i).measurements.Centroid];
                data_temp = data_temp(3:3:end)*cells(i).scaling_dxyz-zOffset;
            otherwise     
                data_temp = [cells(i).measurements.(fieldY)];
        end
        minY = min(data_temp(~isnan(data_temp)));
        maxY = max(data_temp(~isnan(data_temp)));
    end
    rangeY = [minY maxY];
end
binsY = linspace(rangeY(1), rangeY(2), nBinsY);

% Mapping
[X, Y] = meshgrid(binsX, binsY);
Z = cell(size(X));
% Run through data and sort values into bins
for i = 1:numel(cells)
    try
        x = time(i);       
        y = getData(cells(i), fieldY, cells(i).scaling_dxyz, zOffset);
        z = getData(cells(i), fieldZ, cells(i).scaling_dxyz, zOffset);
        
        if islogical(z) 
            z = double(z);
        end
        if islogical(y) 
            y = double(y);
        end
        
        if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(z)
            return;
        end
        
        idxX = find(binsX >= x,1);
        for n = 1:numel(y)
            if y(n) > rangeY(2) || y(n) < rangeY(1) || isnan(y(n))
                continue;
            end
            idxY = find(binsY <= y(n),1, 'last');
            Z{idxY, idxX} = [Z{idxY, idxX} z(n)];
        end
        
    catch err
        fprintf('  - WARNING: "%s"\n', err.message);
    end
end

% Output
try
    eval(['map = cellfun(@(x) ',averagingFcn,'(x(~isnan(x))), Z, ''UniformOutput'', true);']);
catch
    eval(['map = cellfun(@(x) ',averagingFcn,'(x(~isnan(x))), Z, ''UniformOutput'', false);']);
    map = generateUniformOutput(map);
end

emptyEntries = cellfun(@(x) isempty(x), Z, 'UniformOutput', true);
map(emptyEntries) = NaN;
heatmapMatrix = map;

if isempty(rangeZ)
    data_temp = heatmapMatrix(~isnan(heatmapMatrix(:)));
    rangeZ = [min(data_temp) max(data_temp)];
end
        
% Number of data points per tile 
N = cellfun(@(x) numel(x), Z, 'UniformOutput', true);

% Plotting of counts
h = figure('Name', 'Counts');
h_ax = axes('Parent', h);

surf(X, Y, zeros(size(X)), N, 'Parent', h_ax, 'EdgeColor', 'none', 'FaceColor', 'flat')
colormap(h_ax, parula)
view(h_ax, 0,90)
box(h_ax, 'on');
grid(h_ax, 'off');

set(h_ax, 'clim', [min(N(:)) max(N(:))], 'ydir', 'normal', 'ylim', rangeY, 'xlim', rangeX);

xlabel(h_ax, 'time (h)', 'Interpreter', 'none')
ylabel(h_ax, fieldY, 'Interpreter', 'none')

set(h_ax, 'NextPlot', 'add');
c = colorbar;
label = get(c, 'Label');
set(label, 'String', 'counts');
set(label, 'Interpreter', 'none');

% Plotting of data
h = figure('Name', 'Data');
h_ax = axes('Parent', h);

surf(X, Y, zeros(size(X)), heatmapMatrix, 'Parent', h_ax, 'EdgeColor', 'none', 'FaceColor', 'flat')
colormap(h_ax, parula)
view(h_ax, 0,90)
box(h_ax, 'on');
grid(h_ax, 'off');

set(h_ax, 'clim', [rangeZ(1) rangeZ(2)], 'ydir', 'normal', 'ylim', rangeY, 'xlim', rangeX);

xlabel(h_ax, 'time (h)', 'Interpreter', 'none')
ylabel(h_ax, fieldY, 'Interpreter', 'none')

set(h_ax, 'NextPlot', 'add');
c = colorbar;
label = get(c, 'Label');
set(label, 'String', fieldZ);
set(label, 'Interpreter', 'none');
end

%% Additional functions
function output = getData(data, field, scaling, zOffset)
switch field
    case 'z'
        centroids = [data.measurements.Centroid];
        output = centroids(3:3:end)*scaling-zOffset;
    otherwise
        output = [data.measurements.(field)];
end
end

function map = generateUniformOutput(map)
noEntry = cellfun(@(x) isempty(x), map, 'UniformOutput', true);
map(noEntry) = num2cell(repmat(NaN, sum(noEntry(:)),1));
map = cell2mat(map);
end