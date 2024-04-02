function Table = ReadFromFREDOrInputsFolder( SeriesName, DownloadLatestData, Transformation )

    persistent FRED;

    FileName = [ 'Inputs/' SeriesName '.xlsx' ];

    if DownloadLatestData

        if isempty( FRED ) || ~isconnection( FRED )
    
            FRED = fred;
            FRED.DataReturnFormat = 'timetable';
            FRED.DatetimeType = 'datetime';
    
        end
    
        Table = FRED.fetch( SeriesName );
        assert( isscalar( Table.Data ) );
        Table = Table.Data{ 1 };
        assert( size( Table, 2 ) == 1 );
        Table.Properties.VariableNames{ 1 } = SeriesName;
    
        writetimetable( Table, FileName, 'WriteMode', 'replacefile' );

    else

        Table = readtimetable( FileName, 'ReadVariableNames', true, 'VariableNamingRule', 'modify' );

    end

    if ( nargin > 2 ) && ~isempty( Transformation )
        if strcmpi( Transformation, 'Log' )
            Table.( SeriesName ) = log( Table.( SeriesName ) );
            Table.Properties.VariableNames{ 1 } = [ 'Log' SeriesName ];
        elseif strcmpi( Transformation, 'CCRate' )
            Table.( SeriesName ) = 100 * log( 1 + 0.01 * Table.( SeriesName ) );
            Table.Properties.VariableNames{ 1 } = [ 'CCRate' SeriesName ];
        else
            error( 'Unrecognised transformation.' );
        end
    end

end
