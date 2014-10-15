%% Stack Analysis
% This script attempts to identify the location of potential fiber centers
% in a SiC/SiC composite material.
%
% This script requires
%
% * <https://github.com/tonyfast/SpatialStatisticsFFT Spatial Statistics> 
% * <www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave JSON lab>
% * <https://gist.github.com/tonyfast/de3d7b841cea38dc40fa Maximum Gaussian Likelihood estimation>
% * <https://gist.github.com/tonyfast/ffc1c2966f96f98cae6a Fast Radial Feature Detector>
% * <https://gist.github.com/tonyfast/d7f6212f86ee004a4d2b Image Filter Based Peak Finding>
%

%% Add the Codes to the local Path
addpath(genpath('../CleanStats/'))
addpath oldfrst/

stack.name = 'recon_102_2_0_2PIPcure_1p3cm_18keV_ML_1500ms_0to29.tif';
stack.path = '_data';

stack.full = fullfile( stack.path, stack.name );

header = stack;
header.name = horzcat( stack.name, '.json' );
header.full = fullfile( header.path, header.name );

%% Plotting tools

initplot = @(x)close('all');
cleanplot = @(x)set( gcf, 'Position', get(0, 'ScreenSize'));

%% Functions to Normalize Image information

normalize = @(A)( A-min(A(:)) ) ./ ( max(A(:)) - min(A(:)) );
adjust = @(A)reshape(  ... back to original shape
                imadjust( ... adjust image
                        reshape( ... flatten to 2-D image
                                normalize(A), ... normalize from zero to one
                                size(A,1), numel(A)./size(A,1))), ... Reshape to 2-D array
                                size(A) ... Reshape back to original size
                                );

%% Find Centers for all Datasets
% The heavy lifting happens here.
% This script loops over all the available image slices and detects the
% fiber centers.  The values are stored to local mat-files for later.
%
% Once we identify the fiber centers, we will have to center the centers
% from each slice to determine fiber and non-fiber material.

data.image.layers = 30;
layers.id = -2:2;
layers.weight = ones(1, numel(layers.id) );

for ll = ( max(abs(layers.id)) +1 ) : (data.image.layers-max(abs(layers.id)))
    for jj = 1 : numel( layers.id )
        if jj == 1
            A = double( imread( stack.full, ll + layers.id(jj) ) );
        else
            A(:) =  A + layers.weight(jj) .* double( imread( stack.full, ll + layers.id(jj) ) );
        end
    end
%%
% Average the images over multiple layers
    A(:) = A ./ sum( layers.weight );


%% Appromixate the Local of Fibers
% Use Gaussian Likelihood Estimation to product the fiber phases
% Once the fiber phases are found ``poi``, their spatial statistics are
% computed to the approximate the effective radii of the fibers.
%
% The segmentation is used later to identify centers that could be Fibers.
%

%%
% Segmentation
peaks = segmentprob( A, 251,5 );

%%
% The below visualization shows the fiber regions and other phases
% detected.

initplot();
ax(1) = subplot(1,2,1); imshow(normalize(A(:,:,1))');
ax(2) = subplot(1,2,2); pcolor(normalize(peaks.out.phase)');
axis(ax(2),'ij');
axis(ax(2),'equal');
hc = colorbar;
set( get( hc, 'Ylabel'), 'String', 'Possiblity of Circle Center', 'Fontsize',16,'Rotation',270,'VerticalAlignment','Bottom');
colormap jet
linkaxes(ax);
cleanplot();
figure(gcf);

%% Statistics Based Radius Appromixation

initplot();

% Phase of interest from the segmentation
poi = 5; 

cut = 50;
[ F , xx ] = SpatialStatsFFT( peaks.out.phase == poi, [], 'cutoff', cut, ...
                                                          'shift', true);

snapnow;

raddist = [ F( xx.values{1} == 0, : ); ...
        F( :, xx.values{2} == 0 )']';
radmean = mean(raddist,2);
h = plot( -cut: cut, [raddist, radmean], ...
                                        'LineWidth', 3);
xlabel( 'Distance from Fiber Center', ...
                                    'Fontsize', 16);
ylabel( 'Probability of being inside a fiber', ...
                                    'Fontsize', 16);
legend( h, 'x-direction','y-direction','mean');
set( h(end), 'LineStyle',':');
grid on
cleanplot();
figure(gcf);

%% Fast Radial Feature Detector
% The FRFD detectors circular features of a set of prescribed radii to look
% for.  The output image is a transform of the original image.  From the
% transform, centers can be found.

% These radii were chosen using the previous plot
radii = [9 : 12];

T = fastradialv( normalize(A), radii ,2 ); % the third dimension in t corresponds to the input radii
% Average over each radii
T = mean(T,3);


%%
% Illustrate the FRFD

initplot()
ax(1) = subplot(1,2,1); imshow(normalize(A(:,:,1))');
ax(2) = subplot(1,2,2); pcolor(normalize(T)');
axis(ax(2),'ij');
axis(ax(2),'equal');
hc = colorbar;
set( get( hc, 'Ylabel'), 'String', 'Possiblity of Circle Center', 'Fontsize',16,'Rotation',270,'VerticalAlignment','Bottom');
colormap jet
linkaxes(ax); 
cleanplot();
figure(gcf);

%%
close all;
maxim = zeros( size( T ) );
% Find peaks for each radius
for ii = 1 : size( T, 3)
    maxim(:,:,ii) = Find_Peaks( T(:,:,ii) , 'neighborhood', [ 11 11 1], 'diff', true );
end
%% Cluster Peaks from Different Radii

[x,y] = find( maxim );

% ``XX`` is the output
XX  = simplecluster( [x,y], 4 );   
[ xf yf ] =deal( XX(:,1), XX(:,2));

%% Identify center type
I = sub2ind( size( A ), round( XX(:,1) ), round( XX(:,2) ) );
% Precomputed phases
phase = peaks.out.phase( I );
XX(:,3) = phase;

%% Plot the Centers and Their Classes
% Overlay centers and their their class onto the original image
initplot();
co = flipud( cbrewer('qual','Paired',numel(radii)+10 ) );
figure
imshow( normalize(A) );
hold on
for ii = unique( phase )'
    b = phase == ii;
    plot( yf(b), xf(b), 'ko', 'MarkerFaceColor', co(ii,:) );
end

colormap(gray);
hold off
cleanplot();
%% Save information

description = 'Created from stackanalysis.m';
layer = ll;
return
%%
%

% Save the data locally
% save( fullfile('_data','centers', sprintf('clayer_%i.mat', ll ) ), 'XX', 'description', 'layer' );
end