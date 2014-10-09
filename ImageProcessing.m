%% Segmenting Fibers from Cross Sections in CT Data.
% There are several phases in these images and this script solely focuses
% on identifying a good guess for the centers of the fiber positions.
% There are three images from consecutive CT slices.
%
% Requirements: Image Processing Toolbox

%% Import Some Little Widgets
% Normalize and adjust images.  

if ~exist('normalize','var') | ~exist('adjust','var')
    % GIST raw
    rawurl = 'https://gist.githubusercontent.com/tonyfast/8a2bb4752e0cfc55c99f/raw/f706ad03b824c4e17776d012eefd0ec755d133e5/adjust_normalize.m'
    s = urlread(  rawurl );
    eval( s );
    clear( 'rawurl','s')
end

%% For Presentation

setlim = @(x)eval('xlim([60 380]) ;ylim([520 840]);')

%% Import the images

lcldir = '_data';
ims = dir(fullfile( lcldir, '*.tif' ) );

ct = 0;
clear A
for im = ims'
    ct = ct + 1;
    A(:,:,ct) =  double( imresize(...
        imread( fullfile( lcldir, im.name ) ), .5) );
    
end

O = A(:,:,2); 
%% Create Exemplar Image
% The working image is the average of all three slices in the raw images

%% 
% Average All Picture
A = mean(A,3);

%% 
% Normalize Pixel Values
A(:) = normalize(A);

%% 
% Adjust histogram
Aadjust = adjust(A);

%%
% Difference between original and adjusted image.  It is clear that
% different phases have different levels of adjustment in the histogram.
% Could this correspond to different light interactions of the different
% material phases?

dA = Aadjust - A; % Can the difference in pixel adjustment be related to a physical parameters

% Plot %%%%%%%%%%
ax(1) = subplot(1,2,1)
pcolor( A );
axis equal; axis tight; colorbar; shading flat;
title( 'Original Image','Fontsize',16 )
ax(2) = subplot(1,2,2)
pcolor( dA(:,:,1) );
axis equal; axis tight; hc = colorbar; shading flat;
set( get( hc, 'Ylabel'), 'String', 'Normalized Pixel Adjustment' ,...
    'Rotation', 270, 'FontSize', 14,'VerticalAlignment','Bottom')
title( {'Difference After Adjustment', ...
    'Each phase in the image bears a distinct difference after adjustment'} , ...
    'Fontsize',16 )
linkaxes( ax );
figure(gcf)

colormap(cbrewer('div','PuOr',21))

setlim();
%% Identify Fiber Centers 
% Use a recipe of image processing and statistics to find fiber centers and
% classify the different phases the image.

I = normalize( dA ) ;

% Use a Vectorized Fast Radial Transform adopted from:
%
% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% Reference:
% Loy, G.  Zelinsky, A.  Fast radial symmetry for detecting points of
% interest.  IEEE PAMI, Vol. 25, No. 8, August 2003. pp 959-973.


[ I1] = fastradialv( I(:,:,1), 4:6,2);

% Plot %%%%%%%%%%%%%%%%
clf
ax(1) = subplot(1,2,1)
surface( I(:,:,1) ); shading flat; colorbar
axis equal
ax(2) = subplot(1,2,2)
surface( I1(:,:,1) ); shading flat; colorbar
title({'FRST Transform',...
    'High Values indicate circular features'}, ...
    'FontSize',16)
axis equal
linkaxes( ax );
colormap hsv
figure(gcf);
colormap gray
setlim();

%%
% Find the Maxima after the Fast Symmetric Radial Transform has been
% applied.
maxim = Find_Peaks( I1(:,:,1), 'neighborhood', [ 9 9 1] );
[ x, y] = find( maxim );
[ id ] = find( maxim );

clf
ax(1) = subplot(1,2,1)
pcolor( O  ); shading flat; colorbar
hold on
spy( (maxim .* I(:,:,1)) , 10);
hold off
axis equal
ax(2) = subplot(1,2,2)
plot(y,x,'mo','MarkerFaceColor', .1*ones(1,3) )
hold on
[ x2, y2] = deal( y, x);
b = y2 > 500 & y2 < 860 & x2 > 40 & 42 < 400;
voronoi(x2(b),y2(b))
hold off
title('Voronoi Tesselation of Potential Fiber Centers')
axis ij
axis equal
linkaxes( ax );
colormap hsv
figure(gcf);
colormap gray

setlim();
%% Filter out Non-Fiber Phases
% The histogram of the pixel values of the peak centers shows three
% distinct peaks.  We will fit peaks to the histogram then compute the
% probability of each center's phase.

[yy,xx ] = hist(I(id),51);
p = peakfit( yy, 0, 0, 3, 0 );

peaks = struct( 'centerid', p(:,2), ... % mean
    'area', p(:,3), ... 
    'widthid', p(:,4) ); % Standard deviation

figure(gcf);
%%
% Make the histogram peak algorithm a little more structured & normalize
% the signal values.  _I made this harder than it needs to be, whatever_.

peaks = setfield( peaks, 'center', ...
                interp1( 1 : numel(xx), xx, peaks.centerid ) );

peaks = setfield( peaks, 'width', ...
                interp1( 1 : numel(xx), xx, peaks.widthid ) );

%% Classify centers using highest probabilities
% Take the mean and standard deviation from the peak and compute the
% probiblity of each center existing to each peak using a Gaussian.
% The Highest probability wins.

% A value that scales to the proportion
probabilities = exp(bsxfun(@rdivide,...
    -1* (bsxfun( @minus, peaks.center(:)', I(id) ) ).^2,  ...
    (2* (peaks.width(:)').^2 ) ));

[~,idclass] = max( probabilities, [],2);

pcolor( O );
shading flat; axis equal; colorbar;
hold on
co = cbrewer('qual','Dark2',3);
str = {};
for ii = 1 : 3
    b = idclass == ii;
    h(ii) = plot( y(b), x(b), 'o', ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', co(ii,:) );
    str{ii} = sprintf( 'Peak #%i', ii );
end
hold off
figure(gcf)
legend( h, str, 'Fontsize', 16 )
colormap gray
setlim()

centers = struct('x', x, 'y', y, 'id', idclass );

save( fullfile('_data','Test_Difference_Adjusted.mat') , ...
    'centers' );