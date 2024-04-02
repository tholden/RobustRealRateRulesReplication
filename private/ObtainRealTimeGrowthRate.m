function Table = ObtainRealTimeGrowthRate( VariableName, DownloadVintageData )

    % Obtains real time data on the continuously compounded growth rate of the variable VariableName.

    [ HistoricalDates, HistoricalVintages, HistoricalData ] = ObtainHistoricalData( VariableName, DownloadVintageData );

    assert( all( HistoricalDates( 1 : ( end - 1 ) ) + calendarDuration( 0, 1, 0 ) == HistoricalDates( 2 : end ) ) );
    
    assert( issorted( HistoricalVintages ) );

    HistoricalData = 100 * diff( log( HistoricalData ) );
    HistoricalDates( 1 ) = [];
    
    LastObservation = max( isfinite( HistoricalData ) .* ( 1 : size( HistoricalData, 1 ) ).' );

    assert( issorted( LastObservation ) );

    disp( ' ' );

    while true

        MissingVintage = find( LastObservation( 1 : ( end - 1 ) ) + 1 < LastObservation( 2 : end ), 1 );

        if isempty( MissingVintage )
            break
        end

        disp( 'There seems to be a missing vintage between the following two dates:' );
        disp( [ HistoricalVintages( MissingVintage ), HistoricalVintages( MissingVintage + 1 ) ] );
        disp( ' ' );

        % Duplicate the vintage after the missing vintage.
        HistoricalData = [ HistoricalData( :, 1 : MissingVintage ), HistoricalData( :, MissingVintage + 1 ), HistoricalData( :, ( MissingVintage + 1 ) : end ) ];
        HistoricalVintages = [ HistoricalVintages( 1 : MissingVintage ); HistoricalVintages( MissingVintage + 1 ); HistoricalVintages( ( MissingVintage + 1 ) : end ) ];
        LastObservation = [ LastObservation( 1 : MissingVintage ), LastObservation( MissingVintage ) + 1, LastObservation( ( MissingVintage + 1 ) : end ) ];

    end

    CutOffDates = HistoricalDates( LastObservation ) + calendarDuration( 0, 2, 0 );

    PreCutOff = HistoricalVintages <= CutOffDates;

    % Discard vintages which are followed by another with the same last observation date, before the cut-off date.

    Select = [ ( LastObservation( 1 : ( end - 1 ) ) == LastObservation( 2 : end ) ) & PreCutOff( 2 : end ).', false ];

    HistoricalData( :, Select ) = [];
    HistoricalVintages( Select ) = [];
    LastObservation( Select ) = [];

    % Discard vintages which are post cut-off date and which follow another with the same last observation date.

    Select = [ false, LastObservation( 1 : ( end - 1 ) ) == LastObservation( 2 : end ) ];

    HistoricalData( :, Select ) = [];
    HistoricalVintages( Select ) = [];
    LastObservation( Select ) = [];

    assert( all( LastObservation( 1 : ( end - 1 ) ) + 1 == LastObservation( 2 : end ) ) );

    DeleteUpTo = LastObservation( 1 ) - 1;

    HistoricalData( 1 : DeleteUpTo, : ) = [];
    HistoricalDates( 1 : DeleteUpTo, : ) = [];

    assert( size( HistoricalData, 1 ) == size( HistoricalData, 2 ) );
    assert( numel( HistoricalDates ) == numel( HistoricalVintages ) );

    HistoricalData = diag( HistoricalData );

    assert( all( isfinite( HistoricalData ) ) );

    Table = timetable( HistoricalDates, HistoricalVintages, HistoricalData );

    Delays = days( HistoricalVintages - ( HistoricalDates + calendarDuration( 0, 1, 0 ) ) );

    disp( [ 'For variable ' VariableName ' the [median, mean, max] observation delay was:' ] );
    disp( [ median( Delays ), mean( Delays ), max( Delays ) ]  );
    disp( ' ' );

    Delays = days( HistoricalVintages - ( HistoricalDates + calendarDuration( 0, 2, 0 ) ) );

    disp( [ 'For variable ' VariableName ' the number of delays over one month was:' ] );
    disp( sum( Delays > 0 ) );
    disp( ' ' );

end
