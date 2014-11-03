%% Identifying FIBER Centers from its potential candidates
%
% Currently, potential fibers centers have been identified.  The next step
% is to assign when a center is or is not a fiber.
%
% Fibers will be defined as potential fiber centers in contguous layers are
% found nearby to each other.  We assume that in between layers a fiber
% center an not deviate much.
%
% Some potential fiber centers are actually a pyrolytic phase, these can be
% filtered out in the same way.

%% Load in Potential Fiber Centers
%
% From ``stackanalysis.m``, potential fiber centers have been identified
% for the TIF stack
% ``recon_102_2_0_2PIPcure_1p3cm_18keV_ML_1500ms_0to29.tif``.

centers = dir( '_data/centers/*.mat' );

X = [];
for ii = 1 : numel( centers )
    load( fullfile( '_data', 'centers', centers(ii).name ) );
    
    X = [X; horzcat( XX, layer * ones( size(XX,1), 1) )];
end

%%  Remove a Small Section for Testing
% This subvolume tests the efficacy of the current fiber labelling
% approach.  Using a small window it is easier to visualize the results

 w = [ 300 ,500; 600, 700]
b = X(:,1) > w(1,1) & X(:,1) < w(1,2) & ...
    X(:,2) > w(2,1) & X(:,2) < w(2,2) ;
plot3( X(b,1), X(b,2), X(b,4),'k.')
title(' Subvolume Containing Potential Fiber Centers')
axis equal 
axis tight
grid on
set( gcf, 'Position', get(0, 'ScreenSize'))
figure(gcf)

%% Distance Matrix Function
% Computes the distance between "Fiber Centers" in between a reference
% layer and the following layer.
%
% <http://en.wikipedia.org/wiki/Distance_matrix p_2 norm>

dist = @(X,Y)sqrt( ...
    bsxfun( @minus, X(:,1), Y(:,1)' ).^2 +  ...
    bsxfun( @minus, X(:,2), Y(:,2)' ).^2 );


%% Naive Labeling
%
% # Start on First Layer, initialize indices
% # Search for nearest centers in next layer, append prior index.
% # Fibers are indexed from a reference layer to the subsequent layer.

%% Labelling Script

% Euclidean cutoff distance
cut = 4;

% Layers where fiber centers have been found
layers = unique( X(:,end));
% New varaible for the indexed fibers
XX = horzcat( X( b, :) , ...
                zeros( size( X( b, :) , 1), 1 ) );
% Column where the indices go
ncol = size(XX,2);

% Loop over each layer
for ii = [layers'; ....
          1: numel(layers)]    
      
      % Reference Layer
      l1 = XX(:,ncol - 1) == ii(1);
      % Next Layer
      l2 = XX(:,ncol - 1) == ( ii(1) + 1 );
      
      % Points in Reference Layer
      X1 = XX( l1, 1:2);
      % Points in the Next Layer
      X2 = XX( l2, 1:2);
      D = dist( X1, X2);
      
      % Near fiber on next layer relative to the reference layer
      [v,id] = min( D, [],1 );
      
      % Initialize the indices on the first layer
      if ii(2) == 1
          XX( l1 , ncol ) = 1 : sum(l1);
      end
      
      % Previous indices
      prev = XX( l1 , ncol );
      fl2 = find( l2 );
      % Remap prior indices onto the new ones
      XX( fl2( v < cut) , ncol ) = prev(id( v < cut));
      
      % FIrst Pass INitialize indices 
end


%% Plot Indexed Fiber Positions

clear h

% Unique indices 
unind = setdiff( unique( XX(:,ncol) )' ,... ignore zero index
                        0 );

co = rand(100,3);
for ii = [ unind; ...
           1 : numel( unind ) ]
    b = XX(:, ncol ) == ii(1);
    
    C = XX( b, [1 2 4] );
    C = sortrows( C, 3 );
    h(ii) = plot3( C(:,1), C(:,2), C(:,3), 'ko', ...
                                           'MarkerFaceColor', co( ii(2), : ), ...
                                           'MarkerSize', 16 );
    if ii(2) == 1
        hold on;
    end
end
hold off;
grid on

xlabel('Voxels', 'Fontsize',16)
ylabel('Voxels', 'Fontsize',16)
zlabel('Layers', 'Fontsize',16)
daspect( [ 1 1 .25])
title( 'Color Indentifies Continuous Fibers','Fontsize', 16);
set( gcf, 'Position', get(0, 'ScreenSize'))
figure(gcf)
%%