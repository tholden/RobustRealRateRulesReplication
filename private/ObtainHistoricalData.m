function [ HistoricalDates, HistoricalVintages, HistoricalData ] = ObtainHistoricalData( VariableName, DownloadVintageData )

    disp( ' ' );

    if DownloadVintageData

        HistoricalData = FredFetch.vintall( VariableName );
    
        HistoricalDates = datetime( HistoricalData.date, 'ConvertFrom', 'datenum' );
        HistoricalVintages = datetime( HistoricalData.realtime, 'ConvertFrom', 'datenum' );
        HistoricalData = HistoricalData.value;
    
        writematrix( HistoricalDates,    [ 'Inputs/Historical' VariableName '.xlsx' ], 'WriteMode', 'overwritesheet', 'Sheet', 'Dates' );
        writematrix( HistoricalVintages, [ 'Inputs/Historical' VariableName '.xlsx' ], 'WriteMode', 'overwritesheet', 'Sheet', 'Vintages' );
        writematrix( HistoricalData,     [ 'Inputs/Historical' VariableName '.xlsx' ], 'WriteMode', 'overwritesheet', 'Sheet', 'Data' );
    
    else

        HistoricalDates    = readmatrix( [ 'Inputs/Historical' VariableName '.xlsx' ], 'Sheet', 'Dates'   , 'OutputType', 'datetime' );
        HistoricalVintages = readmatrix( [ 'Inputs/Historical' VariableName '.xlsx' ], 'Sheet', 'Vintages', 'OutputType', 'datetime' );
        HistoricalData     = readmatrix( [ 'Inputs/Historical' VariableName '.xlsx' ], 'Sheet', 'Data', 'Range', [ 1, 1, numel( HistoricalDates ), numel( HistoricalVintages ) ] );

    end

end
