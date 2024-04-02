function [ betaGMM, var_betaGMM, t, pt, sigma2, beta2SLS, J, pJ, F, pF ] = hacIV( X, Z, y, beta0 )

    % X first stage variables (potentially including exogenous), Z second stage variables (potentially including exogenous).
    % beta0 is an optional parameter for performing F-tests of beta=beta0. NaNs in beta0 are estimated. If beta0 is not provided, then F is the first stage F-stat.
    % Returns Iterated GMM estimates.

    % Basically following https://www.stata.com/manuals/rivregress.pdf

    n = size( y, 1 );

    assert( size( X, 1 ) == n );
    assert( size( Z, 1 ) == n );

    k = size( X, 2 );
    l = size( Z, 2 );

    [ QZ, ~ ] = qr( Z, 'econ' );

    % ImMZ = QZ * QZ.';

    [ QQZTX, RQZTX ] = qr( QZ.' * X, 'econ' );

    beta2SLS = RQZTX \ ( QQZTX.' * QZ.' * y );

    % Xhat = ( QZ * QZ.' ) * X;

    betaGMM = beta2SLS;

    while true % Iterate the GMM estimator, as better in small samples.

        [ UW, DWinv ] = EstimateCovariance( X, Z, y, betaGMM );
    
        DrootW = diag( 1 ./ ( eps + sqrt( diag( DWinv ) ) ) );
    
        [ QrootWZTX, RrootWZTX ] = qr( DrootW * UW.' * Z.' * X, 'econ' );

        betaGMMOld = betaGMM;
    
        betaGMM = RrootWZTX \ ( QrootWZTX.' * DrootW * UW.' * Z.' * y );

        if all( abs( betaGMM - betaGMMOld ) < 1e-8 )
            break
        end

    end

    [ US, DS ] = EstimateCovariance( X, Z, y, betaGMM );

    DrootS = diag( sqrt( diag( DS ) ) );

    Temp = ( DrootS * US.' * UW * DrootW * QrootWZTX  ) / ( RrootWZTX.' );

    var_betaGMM = ( ( n * n ) / ( n - k ) ) * ( Temp.' * Temp );

    t = betaGMM ./ sqrt( diag( var_betaGMM ) );

    pt = 2 * ( 1 - tcdf( abs( t ), n - k ) );

    u = y - X * betaGMM;

    sigma2 = ( u.' * u ) / ( n - k );

    m = mean( Z .* u );

    DrootInvS = diag( 1 ./ ( eps + diag( DrootS ) ) );

    m = m * US * DrootInvS;

    J = n * ( m * m.' );

    pJ = chi2cdf( J, l - k, 'upper' );

    if nargin > 3

        Select = isfinite( beta0 );

        XA = X( :,  Select );
        XB = X( :, ~Select );

        yB = y - XA * beta0( Select, 1 );

        [ ~, ~, ~, ~, sigma20 ] = hacIV( XB, Z, yB );

        p1 = sum( ~Select );
        p2 = k;

        RSS1 = sigma20 * ( n - p1 );
        RSS2 = sigma2 * ( n - p2 );

        F = ( ( RSS1 - RSS2 ) / ( p2 - p1 ) ) / ( RSS2 / ( n - p2 ) );

        pF = fcdf( F, p2 - p1, n - p2, 'upper' );

    elseif nargout > 8

        beta0 = zeros( l, 1 );
        beta0( std( Z ) < 1e-8 * mean( Z ) ) = NaN;

        F  = zeros( 1, k );
        pF = zeros( 1, k );

        for i = 1 : k
            [ ~, ~, ~, ~, ~, ~, ~, ~, F( k ), pF( k ) ] = hacIV( Z, Z, X( :, k ), beta0 );
        end

    end

end

function [ US, DS ] = EstimateCovariance( X, Z, y, beta )

    u = y - X * beta;

    Zu = Z .* u;

    CZu = Zu( 2 : end, : );
    LZu = Zu( 1 : ( end - 1 ), : );

    SumLZuCZu = sum( LZu .* CZu );
    SumLZuLZu = sum( LZu .* LZu );

    SelectPositive = ( SumLZuLZu > 0 ) | ( ( SumLZuLZu == 0 ) & ( SumLZuCZu >= 0 ) );
    SumLZuLZu( SelectPositive ) = eps + SumLZuLZu( SelectPositive );
    SumLZuLZu( ~SelectPositive ) = SumLZuLZu( ~SelectPositive ) - eps;

    rho = SumLZuCZu ./ SumLZuLZu;

    n = size( Zu, 1 );

    rhoMax = 1 - 1 / sqrt( n ); % Suggested by Sul, Phillips & Choi (2005)

    rho = max( -rhoMax, min( rho, rhoMax ) );

    Zu = CZu - rho .* LZu;

    n = n - 1;

    S = ( 1 / n ) * ( Zu.' * Zu );
    
    q = 2;
    mStar = floor( 20 * ( n / 100 ) ^ ( 4 / 25 ) );
    c_gamma = 2.6614;

    f = Zu;

    SelectMean = std( Z ) < 1e-8 * mean( Z );
    if any( ~SelectMean )
        f( :, SelectMean ) = [];
    end

    f = sum( f, 2 );

    sigma = zeros( mStar + 1, 1 );
    for j = 0 : mStar
        sigma( j + 1 ) = ( 1 / n ) * ( f( ( 1 + j ) : end ).' * f( 1 : ( end - j ) ) );
    end

    sq = 2 * sum( sigma( 2 : end ) .* ( ( 1 : mStar ).' .^ q ) );
    s0 = sigma( 1 ) + 2 * sum( sigma( 2 : end ) );
    gamma = c_gamma * ( ( sq / s0 ) ^ 2 ) ^ ( 1 / ( 2 * q + 1 ) );
    m = gamma * n ^ ( 1 / ( 2 * q + 1 ) );
    m = min( floor( m ), mStar );

    for l = 1 : ( n - 1 )

        CZu = Zu( ( 1 + l ) : end, : );
        LZu = Zu( 1 : ( end - l ), : );

        S = S + ( PZ( l, m ) / n ) * ( CZu.' * LZu + LZu.' * CZu );

    end
    
    S = 0.5 * ( S + S.' );

    Scale = diag( 1 ./ ( 1 - rho ) );

    S = Scale * S * Scale;
    
    S = 0.5 * ( S + S.' );

    [ US, DS ] = schur( S, 'complex' );
    assert( isreal( US ) );
    assert( isreal( DS ) );
    assert( isdiag( DS ) );

end

function K = PZ( l, m )

    z = abs( l ./ ( m + 1 ) );

    K = zeros( size( z ) );
    Select = abs( z ) <= 0.5;
    zSelect = z( Select );
    K( Select ) = 1 - 6 * ( zSelect .* zSelect .* ( 1 - zSelect ) );
    Select = ( ~Select ) & ( abs( z ) <= 1 );
    OMzSelect = 1 - z( Select );
    K( Select ) = 2 * ( OMzSelect .* OMzSelect .* OMzSelect );

end
