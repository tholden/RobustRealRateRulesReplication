function [ A, B, C, D, Mean0, Cov0, StateType ] = PiStarParameterFunction( Parameters )

    % "Hat" below is subscript 1 in the paper!
    % State vector: 1 epsilonHat, 2 epsilonStarHat, 3 PiHat, 4 PiStarHat, 5 PiStarInfinity, 6 Pi0, 7 IMR0LInvPi0, 8 Pi12, 9 IMR12LInvPi12, 10 Pi24, 11 IMR24LInvPi24, 12 Pi36, 13 IMR36LInvPi36, 14 Pi, 15 PiStar, 16-28 LagPiStar1-LagPIStar13
    % Shock vector: 1 epsilonHat, 2 epsilonStarHat, 3 epsilonInfinity, 4 epsilon0, 5 epsilon12, 6 epsilon24, 7 epsilon36
    % Observables vector: 1-40 PiStarAnnual, 41 Pi, 42 PiStarInfinity
    % Parameter vector: 1 Transformed rhoHat, 2 Transformed rhoStarHat, 3 psiHat, 4 psiStarHat, 5 log( sigmaHat ), 6 log( sigmaStarHat ), 7 log( sigmaInfinity ), 8 log( sigma2 )

    rho = @( s ) exp( -1 / ( 1 + s ) );
    scale = @( s ) ( 1 + s ) * rho( s ) ^ s;

    A = diag( [ 0, 0, 2 ./ ( 1 + exp( -Parameters( 1 ) ) ) - 1, 2 ./ ( 1 + exp( -Parameters( 2 ) ) ) - 1, 0.9999, rho( 0 ), rho( 0 ), rho( 12 ), rho( 12 ), rho( 24 ), rho( 24 ), rho( 36 ), rho( 36 ), 0, 0, zeros( 1, 13 ) ] );
    A( 3, 1 ) = Parameters( 3 );
    A( 4, 2 ) = Parameters( 4 );
    A( 15 : 28, 15 : 28 ) = diag( ones( 13, 1 ), -1 );
    
    B = zeros( 28, 7 );
    B( 1, 1 ) = 1;
    B( 2, 2 ) = 1;
    B( 3, 1 ) = exp( Parameters( 5 ) );
    B( 4, 2 ) = exp( Parameters( 6 ) );
    B( 5, 3 ) = exp( Parameters( 7 ) );
    B( 7, 4 ) = exp( Parameters( 8 ) ) / scale( 0 );
    B( 9, 5 ) = exp( Parameters( 8 ) ) / scale( 12 );
    B( 11, 6 ) = exp( Parameters( 8 ) ) / scale( 24 );
    B( 13, 7 ) = exp( Parameters( 8 ) ) / scale( 36 );

    M = eye( 28 ); % Left hand side matrix!
    M( 6, 7 ) = -1;
    M( 8, 9 ) = -1;
    M( 10, 11 ) = -1;
    M( 12, 13 ) = -1;
    M( 14, 3 ) = -1;
    % M( 14, 5 ) = -Parameters( 9 );
    M( 14, 15 ) = -1; % -( 1 - Parameters( 9 ) );
    M( 15, 4 ) = -1;
    M( 15, 5 ) = -1;
    M( 15, 6 ) = -1;
    M( 15, 8 ) = -1;
    M( 15, 10 ) = -1;
    M( 15, 12 ) = -1;

    A = M \ A;
    B = M \ B;

    C = zeros( 42, 28 );
    SelectPiStarAnnual0 = zeros( 1, 28 );
    SelectPiStarAnnual0( 1, 15 : 26 ) = 1;
    SelectPiStarAnnual1 = zeros( 1, 28 );
    SelectPiStarAnnual1( 1, 16 : 27 ) = 1;
    SelectPiStarAnnual2 = zeros( 1, 28 );
    SelectPiStarAnnual2( 1, 17 : 28 ) = 1;
    SelectPiStarAnnual = ( SelectPiStarAnnual0 + SelectPiStarAnnual1 + SelectPiStarAnnual2 ) / 3; % Observed PiStarAnnual is Q4 to Q4 so we need to average over the quarter.
    C( 1, : ) = SelectPiStarAnnual;
    Power = eye( 28 );
    for t = 2 : 40
        Power = A * Power;
        C( t, : ) = SelectPiStarAnnual * Power;
    end
    C( 41, 14 ) = 1;
    C( 42, 5 ) = 1;

    D = zeros( 42, 0 );

    Mean0 = zeros( 28, 1 );
    Mean0( [ 5, 14 : 28 ] ) = 2 / 12;

    Cov0 = [];

    StateType = zeros( 28, 1 );

end
