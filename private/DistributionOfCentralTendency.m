function DistributionOfCentralTendency( nu, Caption, FileName )

    Repeats = 1000;
    
    Mean = zeros( 4, 1 );
    Median = zeros( 4, 1 );
    CentralTendencyMidPoint = zeros( 4, 1 );
    
    parfor Repeat = 1 : Repeats
    
        [ CurrentMean, CurrentMedian, CurrentCentralTendencyMidPoint ] = DistributionOfCentralTendencyWorkChunk( Repeat, nu );
    
        Mean = Mean + CurrentMean;
        Median = Median + CurrentMedian;
        CentralTendencyMidPoint = CentralTendencyMidPoint + CurrentCentralTendencyMidPoint;
    
    end
    
    Moment = ( 1 : 4 ).';
    Mean = Mean / Repeats;
    Median = Median / Repeats;
    CentralTendencyMidPoint = CentralTendencyMidPoint / Repeats;
    
    Table = table( Moment, Mean, Median, CentralTendencyMidPoint );
    
    disp( ' ' );
    disp( Caption );
    disp( ' ' );
    disp( Table );
    disp( ' ' );
    
    writetable( Table, FileName );

end


function [ Mean, Median, CentralTendencyMidPoint ] = DistributionOfCentralTendencyWorkChunk( Seed, nu )

    rng( Seed, 'simdTwister' );
    
    N = 19;
    M = 100000;

    if isfinite( nu )
        T = trnd( nu, N, M );
    else
        T = randn( N, M );
    end

    T = sort( T );
    
    Mean = GetStatistics( mean( T ) );
    Median = GetStatistics( T( 10, : ) );
    CentralTendencyMidPoint = GetStatistics( 0.5 * ( T( 4, : ) + T( 16, : ) ) );

end

function Statistics = GetStatistics( RowVector )
    
    Powers = 1 : 4;
    Statistics = mean( abs( RowVector.' ) .^ Powers ).';

end
