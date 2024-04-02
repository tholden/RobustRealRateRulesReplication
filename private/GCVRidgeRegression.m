function [ YQueryPoints, BestGCVObjective, BetaHat, Lambda, PenaltyScale ] = GCVRidgeRegression( XPoints, YPoints, XQueryPoints, Options )

    Verbosity    = 0;
    Demean       = true;
    Rescale      = true;
    PenaltyScale = 0;
    Debug        = false;

    if nargin == 0
        YQueryPoints = struct( 'Verbosity', Verbosity, 'Demean', Demean, 'Rescale', Rescale, 'PenaltyScale', PenaltyScale, 'Debug', Debug );
        return
    end
    
    if nargin > 3 && isstruct( Options )
        if isfield( Options, 'Verbosity' )
            Verbosity = Options.Verbosity;
        end
        if isfield( Options, 'Demean' )
            Demean = Options.Demean;
        end
        if isfield( Options, 'Rescale' )
            Rescale = Options.Rescale;
        end
        if isfield( Options, 'PenaltyScale' )
            PenaltyScale = Options.PenaltyScale;
        end
        if isfield( Options, 'Debug' )
            Debug = Options.Debug;
        end
    end
    
    Dimension    = size( XPoints, 1 );
    Observations = size( XPoints, 2 );
    QuerySamples = size( XQueryPoints, 2 );
    
    assert( size( XQueryPoints, 1 ) == Dimension, 'GCVRidgeRegression:BadInputNumberOfRows', 'QueryPoints had an unexpected number of rows.' );
    
    assert( numel( YPoints ) == Observations, 'GCVRidgeRegression:DifferentSampleSizes', 'XPoints and YPoints should contain the same number of observations.' );
    
    if isempty( PenaltyScale ) || ~isfinite( PenaltyScale ) || PenaltyScale <= 0
        PenaltyScale = Observations / Dimension;
    end
    
    % Notation follows https://en.wikipedia.org/wiki/Tikhonov_regularization

    X = XPoints.';
    y = YPoints(:);
    
    if Demean
        
        XMean = mean( X );
        X = X - XMean;
        
        yMean = mean( y );
        y = y - yMean;
        
    end
    
    if Rescale
        
        if Demean
            XStd = max( eps, std( X ) );
        else
            XStd = max( eps, sqrt( mean( X .* X ) ) );
        end
        X = X ./ XStd;
        
        if Demean
            yStd = max( eps, std( y ) );
        else
            yStd = max( eps, sqrt( mean( y .* y ) ) );
        end
        y = y ./ yStd;
        
    end
    
    [ U, S, V ] = svd( X, 'econ' );
    
    S = diag( S );
    
    Tolerance = max( size( X ) ) * eps( norm( X ) );
    
    ToDelete = S <= Tolerance;
    
    U( :, ToDelete ) = [];
    S( ToDelete )    = [];
    V( :, ToDelete ) = [];
    
    RSS0 = sum( ( y - U * ( U' * y ) ) .^ 2 );
    
    if isempty( S )
        
        YQueryPoints     = zeros( 1, QuerySamples );
        BestGCVObjective = log( RSS0 ) - 2 * log( Observations );
        BetaHat          = zeros( Dimension, 1 );
        Lambda           = 0;
        
        if Demean
            YQueryPoints = yMean + YQueryPoints;
        end
        
        return
        
    end
    
    [ Lambda, BestGCVObjective ] = GridFMinBound( @( Lambda_ ) FastObjective( Lambda_, U, S, y, RSS0, PenaltyScale ), 0, 1, 21, Verbosity - 1, eps );
    
    if Debug % Paranoid checks that the SVD version really is equivalent to the original.
        
        [ LambdaAlt, BestGCVObjectiveAlt ] = GridFMinBound( @( Lambda_ ) SlowObjective( Lambda_, X, y, PenaltyScale ), 0, 1, 21, Verbosity - 1, eps );
        
        assert( abs( Lambda - LambdaAlt ) < 1e-4 );
        assert( abs( BestGCVObjective - BestGCVObjectiveAlt ) < 1e-4 );
        
        WarningState = warning( 'off', 'MATLAB:fplot:NotVectorized' );
        figure;
        fplot( @( Lambda_ ) SlowObjective( Lambda_, X, y, PenaltyScale ), [ 0.5 * Lambda, 0.5 * ( 1 + Lambda ) ] );
        figure;
        fplot( @( Lambda_ ) FastObjective( Lambda_, U, S, y, RSS0, PenaltyScale ), [ 0.5 * Lambda, 0.5 * ( 1 + Lambda ) ] );
        warning( WarningState );
        
    end
    
    alpha2 = PenaltyScale * Lambda / ( 1 - Lambda );
    
    D = S ./ ( S .* S + alpha2 );
    
    BetaHat = V * diag( D ) * U' * y;
    
    if Debug % Paranoid checks that the SVD version really is equivalent to the original.
        
        BetaHatAlt = inv( X' * X + alpha2 * eye( size( X, 2 ) ) ) * ( X' * y ); %#ok<MINV>
        
        assert( max( abs( BetaHat - BetaHatAlt ) ) < 1e-6 );
        
        if Verbosity > 0
            disp( [ BetaHat BetaHatAlt ] );
        end
        
    end
    
    if Rescale
        BetaHat = BetaHat .* ( YStd ./ XStd.' );
    end
    
    if Demean
        XQueryPoints = [ ones( 1, QuerySamples ); XQueryPoints ];
        BetaHat = [ yMean - XMean * BetaHat; BetaHat ];
    end
    
    YQueryPoints = BetaHat' * XQueryPoints;
    
end

function Objective = SlowObjective( Lambda, X, y, PenaltyScale )

    alpha2 = PenaltyScale * Lambda / ( 1 - Lambda );
    
    [ m, n ] = size( X );
    
    InverseTerm = inv( X' * X + alpha2 * eye( n ) );

    betaHat = InverseTerm * ( X' * y ); %#ok<MINV>
    
    Objective = log( sum( ( X * betaHat - y ) .^ 2 ) ) - 2 * log( trace( eye( m ) - X * InverseTerm * X.' ) ); %#ok<MINV>

end

function Objective = FastObjective( Lambda, U, S, y, RSS0, PenaltyScale )

    alpha2 = PenaltyScale * Lambda / ( 1 - Lambda );
    
    [ m, q ] = size( U );
    
    RatioTerm = alpha2 ./ ( S .* S + alpha2 );

    Objective = log( RSS0 + sum( ( U * ( RatioTerm .* ( U' * y ) ) ) .^ 2 ) ) - 2 * log( m - q + sum( RatioTerm ) );

end
