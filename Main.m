% Replication code for the paper "Robust Real Rate Rules" by Tom D. Holden, to be published in Econometrica.
% Copyright (C) 2024, Tom D. Holden.
% 
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

close all;
clear variables;

DownloadLatestData = false; % Set this to true to update with the latest data, or false to use the saved data from the Inputs folder.

tic;

warning( 'off', 'MATLAB:DELETE:FileNotFound' );
warning( 'off', 'MATLAB:table:ModifiedAndSavedVarnames' );
warning( 'off', 'MATLAB:RMDIR:RemovedFromPath' );
warning( 'off', 'MATLAB:MKDIR:DirectoryExists' );

diary off;
try
    delete Outputs/Log.txt;
catch
end
try
    mkdir Outputs;
catch
end
diary Outputs/Log.txt;

%% Remove Dynare from the path, if present, as it can clash with some Matlab functions (e.g. arima).

DynarePath = which( 'dynare.m' );

if ~isempty( DynarePath )

    % Create a modified version that does not turn on/off the diary!
    Text = readlines( DynarePath );
    Text = regexprep( Text, '^\s*diary\s+\S+\s*;?\s*$', '', 'lineanchors' );
    Text = regexprep( Text, '^\s*diary\s*\(.*\)\s*;?\s*$', '', 'lineanchors' );
    writelines( Text, 'private/dynareNoLog.m' );

    % Get the Dynare root directory.
    DynarePath = fileparts( DynarePath ); % The matlab sub-directory.
    DynarePath = fileparts( DynarePath ); % The root Dynare directory.

    % Remove Dynare from the path.
    PathParts = strsplit( path, ';' );
    Select = ~startsWith( PathParts, DynarePath, 'IgnoreCase', true );
    path( strjoin( PathParts( Select ), ';' ) );

end

%% Create a modified version of autocorr.m.

try
    delete private/autocorrNoMean.m;
catch
end

Path = which( 'autocorr.m' );
Text = readlines( Path );
Text = regexprep( Text, 'y\s*=\s*y\s*-\s*mean\(\s*y.*?\);', '' );
writelines( Text, 'private/autocorrNoMean.m' );

rehash path;

%% Obtain source data, and perform basic initial transformations.

disp( ' ' );
disp( '*** PREPARING DATA ***' );
disp( ' ' );

if DownloadLatestData

    try %#ok<*UNRCH>
        mkdir Inputs;
    catch
    end
    
    websave( 'Inputs/FOMC_Bauer_Swanson.xlsx', 'https://www.michaeldbauer.com/files/FOMC_Bauer_Swanson.xlsx' );

    websave( 'Inputs/SPF.xlsx', 'https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/historical-data/medianlevel.xlsx' );

end

DGS5 = ReadFromFREDOrInputsFolder( 'DGS5', DownloadLatestData );

LogRealPCEPerCapita = ReadFromFREDOrInputsFolder( 'A794RX0Q048SBEA', DownloadLatestData, 'Log' );

LogCPIAUCSL   = ReadFromFREDOrInputsFolder( 'CPIAUCSL', DownloadLatestData, 'Log' );
CCRateT5YIEM  = ReadFromFREDOrInputsFolder(   'T5YIEM', DownloadLatestData, 'CCRate' );
CCRateT5YIFRM = ReadFromFREDOrInputsFolder(  'T5YIFRM', DownloadLatestData, 'CCRate' );

MonthlySourceData = LogCPIAUCSL;
MonthlySourceData = outerjoin( MonthlySourceData, CCRateT5YIEM );
MonthlySourceData = outerjoin( MonthlySourceData, CCRateT5YIFRM );
MonthlySourceData.CPIInflation = [ NaN; 100 * diff( MonthlySourceData.LogCPIAUCSL ) ];

DownloadVintageData = DownloadLatestData && exist( 'api.txt', 'file' );

if DownloadVintageData

    try
        rmdir FredFetch s;
    catch
    end

    gitclone https://github.com/tholden/FredFetch;
    addpath FredFetch;

    copyfile api.txt FredFetch/api.txt f;

elseif DownloadLatestData

    disp( 'To obtain up to date CPI and PCEPI vintages, obtain a FRED API key from https://fredaccount.stlouisfed.org/apikeys and then place it in a file called "api.txt" in this directory.' )

end

RealTimeCCRateCPI   = ObtainRealTimeGrowthRate( 'CPIAUCSL', DownloadVintageData );
RealTimeCCRatePCEPI = ObtainRealTimeGrowthRate(    'PCEPI', DownloadVintageData );

[ HistoricalSEPDates, HistoricalSEPVintages, HistoricalSEPSRData ] = ObtainHistoricalData( 'PCECTPICTM', DownloadVintageData );

CCRatePCECTPICTMLR = ReadFromFREDOrInputsFolder( 'PCECTPICTMLR', DownloadVintageData, 'CCRate' );

HistoricalSEPLRVintages = CCRatePCECTPICTMLR.Properties.RowTimes;
HistoricalSEPLRRawData  = CCRatePCECTPICTMLR.CCRatePCECTPICTMLR;

HistoricalSEPLRData = interp1( HistoricalSEPLRVintages, HistoricalSEPLRRawData, HistoricalSEPVintages, 'nearest', 'extrap' ).'; % Nearest neighbour fill of missing entries (at the start). Also removes one redundant observation in the middle.

HistoricalSEPVintages.Day = 1; % Discard day of the month.

% Now, given CCRatePCECTPICTMLR is mostly constant, we can safely fill in gaps between identical observations to make monthly data. We call this ExtendedHistoricalSEPLRData.

HistoricalSEPLRVintages.Day = 1;

ExtendedHistoricalSEPLRDates = union( RealTimeCCRateCPI.Properties.RowTimes, HistoricalSEPLRVintages );

assert( all( ExtendedHistoricalSEPLRDates( 1 : ( end - 1 ) ) + calendarDuration( 0, 1, 0 ) == ExtendedHistoricalSEPLRDates( 2 : end ) ) );

ExtendedHistoricalSEPLRData = NaN( size( ExtendedHistoricalSEPLRDates ) );

ExtendedHistoricalSEPLRData( ExtendedHistoricalSEPLRDates <= HistoricalSEPLRVintages( 1   ) ) = HistoricalSEPLRRawData( 1   );
ExtendedHistoricalSEPLRData( ExtendedHistoricalSEPLRDates >  HistoricalSEPLRVintages( end ) ) = HistoricalSEPLRRawData( end );

for t = 2 : numel( HistoricalSEPLRRawData )

    if HistoricalSEPLRRawData( t ) == HistoricalSEPLRRawData( t - 1 )
        ExtendedHistoricalSEPLRData( ExtendedHistoricalSEPLRDates <= HistoricalSEPLRVintages( t ) & ExtendedHistoricalSEPLRDates > HistoricalSEPLRVintages( t - 1 ) ) = HistoricalSEPLRRawData( t );
    else
        ExtendedHistoricalSEPLRData( ExtendedHistoricalSEPLRDates == HistoricalSEPLRVintages( t ) ) = HistoricalSEPLRRawData( t );
    end

end

clear HistoricalSEPLRVintages HistoricalSEPLRRawData;

GS5 = ReadFromFREDOrInputsFolder( 'GS5', DownloadLatestData, 'CCRate' );

if DownloadVintageData

    try
        rmdir oilsupplynews s;
    catch
    end

    gitclone https://github.com/dkaenzig/oilsupplynews;

    OilSupplyNewsFiles = dir( 'oilsupplynews' );

    Latest = datetime( 1000, 1, 1 );
    LatestFileName = '';

    for i = 1 : numel( OilSupplyNewsFiles )

        CurrentFileName = OilSupplyNewsFiles( i ).name;

        if endsWith( CurrentFileName, '.xlsx', 'IgnoreCase', true )

            Date = sscanf( CurrentFileName, 'oilSupplyNewsShocks_%dM%d.xlsx' );
            assert( numel( Date ) == 2 );

            Date = datetime( Date( 1 ), Date( 2 ), 1 );

            if Date > Latest
                LatestFileName = CurrentFileName;
                Latest = Date;
            end

        end

    end

    copyfile( [ 'oilsupplynews/' LatestFileName ], 'Inputs/OilSupplyNewsShocksKaenzig.xlsx', 'f' );

    try
        rmdir oilsupplynews s;
    catch
    end

    try
        rmdir FredFetch s;
    catch
    end

end

%% Estimate an ARMA(1,1) model of consumption growth.

disp( ' ' );
disp( '*** ONLINE APPENDIX C ***' );
disp( ' ' );

Dates = LogRealPCEPerCapita.Properties.RowTimes;

disp( 'Sample for the t-distributed consumption growth ARMA(1,1) model (including the observation consumed by differencing) (Y, Q, Y, Q):' );
disp( [ year( Dates( 1 ) ), quarter( Dates( 1 ) ), year( Dates( end ) ), quarter( Dates( end ) ) ] );
disp( ' ' );

disp( 'Estimation of the t-distributed consumption growth ARMA(1,1) model on this sample:' );
disp( ' ' );

Model = arima( 1, 0, 1 );
Model.Distribution = 't';

estimate( Model, diff( LogRealPCEPerCapita.LogA794RX0Q048SBEA ), 'AR0', 0.6, 'MA0', -0.4, 'Display', 'full', 'Options', optimoptions( 'fmincon', 'ConstraintTolerance', 1e-12, 'Display', 'iter-detailed', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 1e12, 'MaxIterations', 1e12, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12, 'UseParallel', true ) );

Select = year( Dates ) <= 2019;

LogRealPCEPerCapita = LogRealPCEPerCapita( Select, : );

Dates = LogRealPCEPerCapita.Properties.RowTimes;

disp( 'Sample for the Gaussian-distributed consumption growth ARMA(1,1) model (including the observation consumed by differencing) (Y, Q, Y, Q):' );
disp( [ year( Dates( 1 ) ), quarter( Dates( 1 ) ), year( Dates( end ) ), quarter( Dates( end ) ) ] );
disp( ' ' );

disp( 'Estimation of the Gaussian-distributed consumption growth ARMA(1,1) model on this sample:' );
disp( ' ' );

Model = arima( 1, 0, 1 );

estimate( Model, diff( LogRealPCEPerCapita.LogA794RX0Q048SBEA ), 'AR0', 0.6, 'MA0', -0.4, 'Display', 'full', 'Options', optimoptions( 'fmincon', 'ConstraintTolerance', 1e-12, 'Display', 'iter-detailed', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 1e12, 'MaxIterations', 1e12, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12, 'UseParallel', true ) );

%% Examine how close (in some sense) the monetary rule of Smets Wouters (2007) is to a real rate rule.

if ~isempty( DynarePath )

    disp( 'Examination of how close (in some sense) the monetary rule of Smets Wouters (2007) is to a real rate rule.' );
    disp( ' ' );

    addpath( [ DynarePath '/matlab' ] );

    cd SmetsWouters2007;

    dynareNoLog Smets_Wouters_2007_45.mod noclearall nolog console nograph nointeractive nopreprocessoroutput nowarn notime;

    cd ..;

    disp( ' ' );
    disp( 'Covariance matrix of z and the real interest rate:' );
    disp( oo_.var );
    disp( ' ' );
    disp( 'Correlation between z and the real interest rate:' );
    disp( oo_.var( 1, 2 ) / sqrt( oo_.var( 1, 1 ) * oo_.var( 2, 2 ) ) );
    disp( ' ' );
    disp( 'Standard deviation of z and the real interest rate:' );
    disp( sqrt( diag( oo_.var ) ).' );
    disp( ' ' );

    % Remove Dynare from the path again.
    PathParts = strsplit( path, ';' );
    Select = ~startsWith( PathParts, DynarePath, 'IgnoreCase', true );
    path( strjoin( PathParts( Select ), ';' ) );

end

%% Do SPF forecasts predict breakeven inflation?

disp( ' ' );
disp( '*** SUPPLEMENTAL APPENDIX I.1 ***' );
disp( ' ' );

CPI5YR = readmatrix( 'Inputs/SPF.xlsx', 'Sheet', 'CPI5YR' );
CPIF5  = readmatrix( 'Inputs/SPF.xlsx', 'Sheet',  'CPIF5' );

assert( size( CPI5YR, 1 ) == size( CPIF5, 1 ) );
assert( all( CPI5YR( :, 1 : 2 ) == CPIF5( :, 1 : 2 ), 'all' ) );

SPFYears = CPI5YR( :, 1 );
SPFQuarters = CPI5YR( :, 2 );

SPFForecast = 100 * log( 1 + 0.01 * CPI5YR( :, 3 ) );
SPFForwardForecast = 100 * log( 1 + 0.01 * CPIF5( :, 3 ) );

Select = isfinite( SPFForecast ) | isfinite( SPFForwardForecast );

SPFYears    = SPFYears( Select );
SPFQuarters = SPFQuarters( Select );
SPFForecast = SPFForecast( Select );
SPFForwardForecast = SPFForwardForecast( Select );

disp( 'Sample (Y, Q, Y, Q):' )
disp( [ SPFYears( 1 ), SPFQuarters( 1 ), SPFYears( end ), SPFQuarters( end ) ] );

LastMonthKnownBySPFParticipants = datetime( SPFYears, ( SPFQuarters - 1 ) * 3, 1 );
CorrespondingBreakevenDate = LastMonthKnownBySPFParticipants + calendarDuration( 0, 1, 0 );
FirstMonthCoveredInBreakeven = CorrespondingBreakevenDate - calendarDuration( 0, 3, 0 );
FirstMonthCoveredInForecast = datetime( SPFYears - 1, 11, 1 );
FirstMonthCoveredInModifiedForecast = max( FirstMonthCoveredInBreakeven, FirstMonthCoveredInForecast );
assert( all( ( FirstMonthCoveredInModifiedForecast == FirstMonthCoveredInBreakeven ) | ( FirstMonthCoveredInModifiedForecast == FirstMonthCoveredInBreakeven + calendarDuration( 0, 1, 0 ) ) ) );

ModifiedBreakevenInflation = MonthlySourceData.CCRateT5YIEM( CorrespondingBreakevenDate );
Select = FirstMonthCoveredInModifiedForecast > FirstMonthCoveredInBreakeven;
assert( all( Select == ( SPFQuarters == 1 ) ) );
ModifiedBreakevenInflation( Select ) = ModifiedBreakevenInflation( Select ) + ( 6 * MonthlySourceData.CCRateT5YIFRM( CorrespondingBreakevenDate( Select ) ) - 3 * MonthlySourceData.CPIInflation( FirstMonthCoveredInBreakeven( Select ) ) - 2 * MonthlySourceData.CPIInflation( FirstMonthCoveredInBreakeven( Select ) + calendarDuration( 0, 1, 0 ) ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInBreakeven( Select ) + calendarDuration( 0, 2, 0 ) ) ) * ( 1 / 180 );

ModifiedSPFForecast = SPFForecast;
Select = SPFQuarters >= 2;
ModifiedSPFForecast( Select ) = ModifiedSPFForecast( Select ) + ( SPFForwardForecast( Select ) - ( 1 / 3 ) * MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) ) - ( 2 / 3 ) * MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 1, 0 ) ) ) * ( 1 / 60 );
Select = SPFQuarters >= 3;
ModifiedSPFForecast( Select ) = ModifiedSPFForecast( Select ) + ( 3 * SPFForwardForecast( Select ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 2, 0 ) ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 3, 0 ) ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 4, 0 ) ) ) * ( 1 / 60 );
Select = SPFQuarters == 4;
ModifiedSPFForecast( Select ) = ModifiedSPFForecast( Select ) + ( 3 * SPFForwardForecast( Select ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 5, 0 ) ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 6, 0 ) ) - MonthlySourceData.CPIInflation( FirstMonthCoveredInForecast( Select ) + calendarDuration( 0, 7, 0 ) ) ) * ( 1 / 60 );

y = ModifiedBreakevenInflation;
x = ModifiedSPFForecast;

[ theta, var_theta, ~, Ptheta0, ~, ~, ~, ~, ~, Ptheta1 ] = hacIV( [ ones( size( x ) ), x ], [ ones( size( x ) ), x ], y, [ NaN; 1 ] );

disp( 'Results of predicting five-year breakeven inflation from five-year SPF inflation expectations (slope, standard error, p-value of == 0, p-value of == 1):' );
disp( [ theta( 2 ) sqrt( var_theta( 2, 2 ) ) Ptheta0( 2 ), Ptheta1 ] );

figure;
scatter( x, y, 12, 'filled', 'k' );
hold on;
XLim = xlim;
YLim = ylim;
plot( XLim, XLim * theta( 2 ) + theta( 1 ), 'k' );
xlabel( 'Five-year SPF inflation expectations (annual rate)' );
ylabel( 'Five-year breakeven inflation (annual rate)')
xlim( XLim );
ylim( YLim );
set( gca, 'Box', 'on' );
saveas( gcf, 'Outputs/SupplementalFigure1BreakevenInflationVsSPF.svg' );

%% Does breakeven inflation forecast actual inflation?

disp( ' ' );
disp( '*** SUPPLEMENTAL APPENDIX I.2 ***' );
disp( ' ' );

PaddedLogCPIAUCSL = [ MonthlySourceData.LogCPIAUCSL; NaN( 60, 1 ) ];

Realised5YearInflation = [ NaN; 20 * ( PaddedLogCPIAUCSL( 61 : end ) - PaddedLogCPIAUCSL( 1 : ( end - 60 ) ) ) ]; % 20 = 100 / 5.
Realised5YearInflation = Realised5YearInflation( 1 : size( MonthlySourceData, 1 ) );

MonthlySourceData.Realised5YearInflation = Realised5YearInflation;

BreakevenAtEffectiveDate = CCRateT5YIEM;

BreakevenAtEffectiveDate.Time = BreakevenAtEffectiveDate.Time - calendarDuration( 0, 3, 0 );
BreakevenAtEffectiveDate.Properties.VariableNames{ 1 } = 'BreakevenAtEffectiveDate';

MonthlySourceData = outerjoin( MonthlySourceData, BreakevenAtEffectiveDate );

Select = isfinite( MonthlySourceData.Realised5YearInflation ) & isfinite( MonthlySourceData.BreakevenAtEffectiveDate );
Times = MonthlySourceData.Time( Select );

disp( 'Sample:' )
disp( [ Times( 1 ), Times( end ) ]  + calendarDuration( 0, 3, 0 ) );

x = MonthlySourceData.BreakevenAtEffectiveDate( Select );
y = MonthlySourceData.Realised5YearInflation( Select );

[ theta, var_theta, ~, Ptheta ] = hacIV( [ ones( size( x ) ), x ], [ ones( size( x ) ), x ], y );

disp( 'Results of predicting five-year realised inflation from five-year breakeven inflation (slope, standard error, p-value):' );
disp( [ theta( 2 ) sqrt( var_theta( 2, 2 ) ) Ptheta( 2 ) ] );

x = diff( x );
y = diff( y );

[ theta, var_theta, ~, Ptheta ] = hacIV( x, x, y );

disp( 'Results of predicting changes in five-year realised inflation from changes in five-year breakeven inflation (slope, standard error, p-value):' );
disp( [ theta sqrt( var_theta ) Ptheta ] );

figure;
scatter( x, y, 12, 'filled', 'k' );
hold on;
XLim = xlim;
YLim = ylim;
plot( XLim, XLim * theta, 'k' );
xlabel( 'Change in five-year breakeven inflation (annual rate)' );
ylabel( 'Change in five-year realised inflation (annual rate)')
xlim( XLim );
ylim( YLim );
set( gca, 'Box', 'on' );
saveas( gcf, 'Outputs/SupplementalFigure2BreakevenInflationVsReality.svg' );

%% Properties of the mid-point of the central tendency.

disp( ' ' );
disp( '*** SUPPLEMENTAL APPENDIX J.1 ***' );
disp( ' ' );

DistributionOfCentralTendency(   5, 'Supplemental Table 1:', 'Outputs/SupplementalTable1DistributionOfCentralTendencyT.xlsx' );
DistributionOfCentralTendency( Inf, 'Supplemental Table 2:', 'Outputs/SupplementalTable2DistributionOfCentralTendencyN.xlsx' );

%% PCE to CPI conversion.

disp( ' ' );
disp( '*** SUPPLEMENTAL APPENDIX J.2 ***' );
disp( ' ' );

gcp;

Data = innerjoin( RealTimeCCRateCPI, RealTimeCCRatePCEPI );

% Remove unneeded pre-2002 data.
Data( year( Data.Properties.RowTimes ) < 2002, : ) = [];

Dates = Data.Properties.RowTimes;

T = size( Data, 1 );

assert( all( Dates( 1 : ( end - 1 ) ) + calendarDuration( 0, 1, 0 ) == Dates( 2 : end ) ) );

disp( 'Sample for PCE to CPI conversion:' );
disp( [ Dates( 1 ), Dates( end ) ] );
disp( ' ' );

Delays = days( Data.HistoricalVintages_RealTimeCCRateCPI - ( Dates + calendarDuration( 0, 1, 0 ) ) );

disp( 'In the merged sample, for CPI the [median, mean, max] observation delay was:' );
disp( [ median( Delays ), mean( Delays ), max( Delays ) ]  );
disp( ' ' );

Delays = days( Data.HistoricalVintages_RealTimeCCRateCPI - ( Dates + calendarDuration( 0, 2, 0 ) ) );

disp( 'In the merged sample, for CPI the number of delays over one month was (compared to the size of the sample):' );
disp( [ sum( Delays > 0 ), T ] );
disp( ' ' );

Delays = days( Data.HistoricalVintages_RealTimeCCRatePCEPI - ( Dates + calendarDuration( 0, 1, 0 ) ) );

disp( 'In the merged sample, for PCEPI the [median, mean, max] observation delay was:' );
disp( [ median( Delays ), mean( Delays ), max( Delays ) ]  );
disp( ' ' );

Delays = days( Data.HistoricalVintages_RealTimeCCRatePCEPI - ( Dates + calendarDuration( 0, 2, 0 ) ) );

disp( 'In the merged sample, for PCEPI the number of delays over one month was (compared to the size of the sample):' );
disp( [ sum( Delays > 0 ), T ] );
disp( ' ' );

CPI = Data.HistoricalData_RealTimeCCRateCPI;
PCEPI = Data.HistoricalData_RealTimeCCRatePCEPI;

% Get initial values using a ridge regression approach.

X0 = [ ones( size( PCEPI ) ), PCEPI ];

alpha_beta0 = regress( CPI, X0 );

CPI0 = CPI - X0 * alpha_beta0;

Triangle = tril( ones( T ) );

Triangle( :, 1 ) = 1000 * T * Triangle( :, 1 ); % The 1000 * T * effectively removes the regularization of the constant.

XPoints = [ Triangle, 6 * diag( PCEPI ) * Triangle ].'; % 6 = 12 (months in a year) / 2 (the inflation target)

[ YQueryPoints, ~, BetaHat ] = GCVRidgeRegression( XPoints, CPI0.', XPoints, struct( 'Demean', false, 'Rescale', false, 'Debug', false, 'Verbosity', 10, 'PenaltyScale', 200 ) );

alpha = alpha_beta0( 1 ) + 1000 * T * BetaHat( 1 ) + cumsum( [ 0; BetaHat( 2 : T ) ] );
beta = alpha_beta0( 2 ) + 6000 * T * BetaHat( T + 1 ) + cumsum( [ 0; 6 * BetaHat( ( T + 2 ) : end ) ] );

Predicted = alpha + beta .* PCEPI;
AltPredicted = YQueryPoints.' + X0 * alpha_beta0;

disp( 'Max error from alternative representation:' );
disp( max( abs( Predicted - AltPredicted ) ) );
disp( ' ' );

Residuals = CPI - Predicted;

Model = ssm( eye( 2 ), [ NaN, NaN; 0, NaN ], arrayfun( @( t ) [ 1, PCEPI( t ) ], 1 : T, 'UniformOutput', false ), NaN, 'Mean0', [ 0; 2 / 12 ], 'Cov0', diag( [ 1, 12 / 2 ] .^ 2 ), 'StateType', [ 2; 2 ] );

CholCov = chol( cov( [ diff( beta ), diff( alpha ) ] ), 'lower' );

OriginalParameters = [ CholCov( 2, 2 ); CholCov( 2, 1 ); CholCov( 1, 1 ); std( Residuals ) ];

[ EstimatedModel, EstimatedParameters, EstimatedParametersCovariance ] = estimate( Model, CPI, OriginalParameters, 'CovMethod', 'sandwich', 'Display', 'full', 'SquareRoot', true, 'Options', optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter-detailed', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 1e12, 'MaxIterations', 1e12, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', true ) );

disp( ' ' );
disp( 'Initial and final parameters:' )
disp( [ OriginalParameters, EstimatedParameters ] );
disp( ' ' );


disp( 'Final parameters in natural domain, with their standard error:' );
ParamterNames = { 'sigma_alpha_alpha'; 'sigma_alpha_beta'; 'sigma_beta_beta'; 'sigma_epsilon' };
StandardErrors = sqrt( diag( EstimatedParametersCovariance ) );
disp( table( EstimatedParameters, StandardErrors, 'RowNames', ParamterNames ) );
disp( ' ' );


alpha_beta = smooth( EstimatedModel, CPI, 'SquareRoot', true );

alpha = alpha_beta( :, 1 );
beta = alpha_beta( :, 2 );

Predicted = alpha + beta .* PCEPI;

stdPCE = std( PCEPI );
beta_stdPCE = beta * stdPCE;
rho_PCE_CPI = beta_stdPCE ./ sqrt( beta_stdPCE .* beta_stdPCE + EstimatedParameters( 3 ) ^ 2 );

Residuals = CPI - Predicted;

figure;
plot( Dates, CPI, Dates, Predicted );
figure;
plot( Dates, Residuals );
figure;
plot( Dates, alpha, 'k' );
xlabel( 'Year' );
ylabel( '$\alpha_t$', 'Interpreter', 'latex' );
saveas( gcf, 'Outputs/SupplementalFigure3alpha.svg' );

figure;
plot( Dates, beta, 'k' );
xlabel( 'Year' );
ylabel( '$\beta_t$', 'Interpreter', 'latex' );
saveas( gcf, 'Outputs/SupplementalFigure4beta.svg' );

figure;
plot( Dates, rho_PCE_CPI, 'k' );
xlabel( 'Year' );
ylabel( 'Implied correlation between PCEPI and CPI inflation' );
saveas( gcf, 'Outputs/SupplementalFigure5rho_PCEPI_CPI.svg' );

                % writetimetable( timetable( Dates, alpha, beta ), 'Outputs/EstimatedAlphaBeta.xlsx', 'WriteMode', 'replacefile' );

% Convert ExtendedHistoricalSEPLRData, HistoricalSEPLRData and HistoricalSEPSRData to CPI

AlphaBetaDates = Dates + calendarDuration( 0, 1, 0 ); % AlphaBetaDates has the month at which this is observed.

alphaAnnual = 12 * alpha;

for t = 1 : numel( ExtendedHistoricalSEPLRDates )

    CurrentDate = ExtendedHistoricalSEPLRDates( t );

    Index = find( AlphaBetaDates <= CurrentDate, 1, 'last' );

    if isscalar( Index )
        ExtendedHistoricalSEPLRData( t ) = alphaAnnual( Index ) + beta( Index ) * ExtendedHistoricalSEPLRData( t );
    else
        ExtendedHistoricalSEPLRData( t ) = NaN;
    end

end

T = numel( HistoricalSEPVintages );

for t = 1 : T

    CurrentDate = HistoricalSEPVintages( t );

    Index = find( AlphaBetaDates <= CurrentDate, 1, 'last' );
    assert( isscalar( Index ) );

    HistoricalSEPLRData(    t ) = alphaAnnual( Index ) + beta( Index ) * HistoricalSEPLRData(    t );
    HistoricalSEPSRData( :, t ) = alphaAnnual( Index ) + beta( Index ) * HistoricalSEPSRData( :, t );

end

%% Estimating the monetary rule on annual data.

disp( ' ' );
disp( '*** SUPPLEMENTAL APPENDIX J.3 ***' );
disp( ' ' );

% Estimate persistence of forecasts.

disp( 'Sample for forecast persistence estimation:' );
disp( HistoricalSEPVintages( [ 1, end ] ).' );
disp( ' ' );

RelativeHistoricalSEPSRData = HistoricalSEPSRData - HistoricalSEPLRData;

x = NaN( T, 1 );
y = NaN( T, 1 );

for t = 1 : T

    Index = find( isfinite( RelativeHistoricalSEPSRData( :, t ) ), 1, 'last' );
    assert( ~isempty( Index ) );
    assert( isfinite( RelativeHistoricalSEPSRData( Index - 1, t ) ) );

    x( t ) = RelativeHistoricalSEPSRData( Index - 1, t );
    y( t ) = RelativeHistoricalSEPSRData( Index    , t );

end

[ rhoRelativeInflation, var_rhoRelativeInflation ] = hacIV( x, x, y );

disp( 'Estimated persistence (and standard error):' )
disp( [ rhoRelativeInflation, sqrt( var_rhoRelativeInflation ) ] );
disp( ' ' );

% Project forward using the simple AR(1) model.

ExtendedHistoricalSEPSRData = [ HistoricalSEPSRData; NaN( 3, T ) ];

for t = 2 : size( ExtendedHistoricalSEPSRData, 1 )

    Select = isfinite( ExtendedHistoricalSEPSRData( t - 1, : ) ) & ~isfinite( ExtendedHistoricalSEPSRData( t, : ) );
    ExtendedHistoricalSEPSRData( t, Select ) = HistoricalSEPLRData( Select ) + rhoRelativeInflation * ( ExtendedHistoricalSEPSRData( t - 1, Select ) - HistoricalSEPLRData( Select ) );

end

% Aggregate forecasts to annual.

YearHistoricalSEPVintages = year( HistoricalSEPVintages );
YearSetHistoricalSEPVintages = unique( YearHistoricalSEPVintages );

CompleteYears = YearSetHistoricalSEPVintages( histcounts( categorical( YearHistoricalSEPVintages, YearSetHistoricalSEPVintages ) ) >= 3 );

disp( 'Sample for annual estimation of RRR:' );
disp( [ CompleteYears( 1 ) + 1, CompleteYears( end ) ] );
disp( ' ' );

AnnualHistoricalSEPSRData = NaN( size( ExtendedHistoricalSEPSRData, 1 ), numel( CompleteYears ) );
AnnualHistoricalSEPLRData = NaN( 1, numel( CompleteYears ) );

for t = 1 : numel( CompleteYears )

    Year = CompleteYears( t );

    Indices = find( YearHistoricalSEPVintages == Year );
    assert( numel( Indices ) >= 3 );

    AnnualHistoricalSEPSRData( :, t ) = mean( ExtendedHistoricalSEPSRData( :, Indices ), 2 );
    AnnualHistoricalSEPLRData( t ) = mean( HistoricalSEPLRData( Indices ) );

end

% Construct current and projections of PiStar

PiStarCurrent    = NaN( numel( CompleteYears ), 1 );
PiStarProjection = NaN( numel( CompleteYears ), 1 );

YearHistoricalSEPSRDates = year( HistoricalSEPDates );

for t = 1 : numel( CompleteYears )

    Year = CompleteYears( t );

    Index = find( YearHistoricalSEPSRDates == Year );
    assert( isscalar( Index ) );
    
    PiStarCurrent( t ) = AnnualHistoricalSEPSRData( Index, t );

    PiStarProjection( t ) = mean( AnnualHistoricalSEPSRData( Index + ( 1 : 5 ), t ) );

end

DecemberCPIIndices = find( month( LogCPIAUCSL.Properties.RowTimes ) == 12 & year( LogCPIAUCSL.Properties.RowTimes ) >= CompleteYears( 1 ) & year( LogCPIAUCSL.Properties.RowTimes ) <= CompleteYears( end ) );

assert( all( month( LogCPIAUCSL.Properties.RowTimes( DecemberCPIIndices + 1  ) ) == 1  ) );
assert( all( month( LogCPIAUCSL.Properties.RowTimes( DecemberCPIIndices - 12 ) ) == 12 ) );
assert( all( month( LogCPIAUCSL.Properties.RowTimes( DecemberCPIIndices - 11 ) ) == 1  ) );

assert( numel( LogCPIAUCSL.Properties.RowTimes( DecemberCPIIndices ) ) == numel( CompleteYears ) );
assert( all( year( LogCPIAUCSL.Properties.RowTimes( DecemberCPIIndices ) ) == CompleteYears ) );

AnnualCPIInflationData = 50 * ( LogCPIAUCSL.LogCPIAUCSL( DecemberCPIIndices ) - LogCPIAUCSL.LogCPIAUCSL( DecemberCPIIndices - 12 ) + LogCPIAUCSL.LogCPIAUCSL( DecemberCPIIndices + 1 ) - LogCPIAUCSL.LogCPIAUCSL( DecemberCPIIndices - 11 ) );

DecemberBreakevenIndices = find( month( CCRateT5YIEM.Properties.RowTimes ) == 12 & year( CCRateT5YIEM.Properties.RowTimes ) >= CompleteYears( 1 ) & year( CCRateT5YIEM.Properties.RowTimes ) <= CompleteYears( end ) );

assert( all( month( CCRateT5YIEM.Properties.RowTimes( DecemberBreakevenIndices + 1 ) ) == 1 ) );

assert( numel( CCRateT5YIEM.Properties.RowTimes( DecemberBreakevenIndices ) ) == numel( CompleteYears ) );
assert( all( year( CCRateT5YIEM.Properties.RowTimes( DecemberBreakevenIndices ) ) == CompleteYears ) );

AnnualBreakevenData = 0.5 * ( CCRateT5YIEM.CCRateT5YIEM( DecemberBreakevenIndices ) + CCRateT5YIEM.CCRateT5YIEM( DecemberBreakevenIndices + 1 ) );

RelativeBreakevenInflation = AnnualBreakevenData - PiStarProjection;

x = AnnualCPIInflationData( 2 : end ) - PiStarCurrent( 2 : end );
y = RelativeBreakevenInflation( 2 : end ) - RelativeBreakevenInflation( 1 : ( end - 1 ) );
Labels = CompleteYears( 2 : end );

[ theta, ~, ~, Ptheta ] = hacIV( x, x, y );

disp( 'Estimated theta from the simple annual exercise:' );
disp( theta );
disp( 'HAC P-value from the simple annual exercise:' );
disp( Ptheta );
disp( 'R^2 from the simple annual exercise:' );
disp( 1 - sum( ( y - x * theta ) .^ 2 ) ./ sum( y .* y ) );
disp( ' ' );

figure;
scatter( x, y, 'filled', 'k' );
hold on;
OnLeft = ismember( Labels, [ 2020, 2011, 2018, 2021 ] );
OnRightSpecial = ismember( Labels, 2010 );
OnRight = ~( OnLeft | OnRightSpecial );
text( x( OnRightSpecial ) + 0.05, y( OnRightSpecial ) + 0.05, num2str( Labels( OnRightSpecial ) ), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( x( OnRight ) + 0.05, y( OnRight ), num2str( Labels( OnRight ) ), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( x( OnLeft ) - 0.05, y( OnLeft ), num2str( Labels( OnLeft ) ), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle' );
XLim = xlim;
YLim = ylim;
plot( XLim, XLim * theta, 'k' );
xlim( XLim );
ylim( YLim );
set( gca, 'Box', 'on' );
saveas( gcf, 'Outputs/SupplementalFigure6AnnualExercise.svg' );

%% Estimating a monthly model for PiStar.

disp( ' ' );
disp( '*** SUPPLEMENTAL APPENDIX J.4 ***' );
disp( ' ' );

Select = year( ExtendedHistoricalSEPLRDates ) >= 2007;

Dates = ExtendedHistoricalSEPLRDates( Select );
PiStarInfinityObs = ExtendedHistoricalSEPLRData( Select ) / 12;

T = numel( Dates );

disp( 'Sample for PiStar estimation:' );
disp( Dates( [ 1, end ] ).' );
disp( ' ' );

PiObs = NaN( T, 1 );

PiObs( ismember( Dates, RealTimeCCRateCPI.Properties.RowTimes ) ) = RealTimeCCRateCPI.HistoricalData( ismember( RealTimeCCRateCPI.Properties.RowTimes, Dates ) );

PiStarAnnualObs = NaN( T, 40 );

for t = 1 : numel( HistoricalSEPVintages )

    VintageDate = HistoricalSEPVintages( t );

    CurrentDates = VintageDate - calendarDuration( 0, 12, 0 ) + calendarDuration( 0, 1 : 40, 0 );

    IsFiniteHistoricalSEPSRData = isfinite( HistoricalSEPSRData( :, t ) );

    FoundDates = HistoricalSEPDates( IsFiniteHistoricalSEPSRData );

    [ IsMember, IndexIntoCurrentDates ] = ismember( FoundDates, CurrentDates );

    assert( all( IsMember ) );

    Index = find( Dates == VintageDate );
    assert( isscalar( Index ) );

    PiStarAnnualObs( Index, IndexIntoCurrentDates ) = HistoricalSEPSRData( IsFiniteHistoricalSEPSRData, t ).';

end

% Obtain initial values from an ARMA(1,1) model for PiObs.

Model = arima( 1, 0, 1 );

EstimatedModel = estimate( Model, PiObs, 'AR0', 0.9, 'MA0', -0.8, 'Display', 'full', 'Options', optimoptions( 'fmincon', 'ConstraintTolerance', 1e-12, 'Display', 'iter-detailed', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 1e12, 'MaxIterations', 1e12, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12, 'UseParallel', true ) );

AR = EstimatedModel.AR{ 1 };
MA = EstimatedModel.MA{ 1 };
Sigma = sqrt( EstimatedModel.Variance );

LogitAR = log( 1 + AR ) - log( 1 - AR );

% Estimate the full model for PiStar.

OriginalParameters = [ LogitAR; LogitAR; MA * Sigma / 2; MA * Sigma / 2; log( Sigma / 2 ); log( Sigma / 2 ); -5; -5 ];

Model = ssm( @PiStarParameterFunction );

Data = [ PiStarAnnualObs, PiObs, PiStarInfinityObs ];

[ EstimatedModel, EstimatedParameters, EstimatedParametersCovariance ] = estimate( Model, Data, OriginalParameters, 'CovMethod', 'sandwich', 'Display', 'full', 'SquareRoot', true, 'Options', optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter-detailed', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', 1e12, 'MaxIterations', 1e12, 'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'UseParallel', true ) );

disp( ' ' );
disp( 'Initial and final parameters:' )
disp( [ OriginalParameters, EstimatedParameters ] );
disp( ' ' );

NaturalDomainParameters = [ 2 ./ ( 1 + exp( -EstimatedParameters( 1 : 2 ) ) ) - 1; EstimatedParameters( 3 : 4 ); exp( EstimatedParameters( 5 : 8 ) ) ];
DNaturalDomainParameters = [ 2 .* exp( -EstimatedParameters( 1 : 2 ) ) ./ ( 1 + exp( -EstimatedParameters( 1 : 2 ) ) ) .^ 2; ones( 2, 1 ); exp( EstimatedParameters( 5 : 8 ) ) ];
disp( 'Final parameters in natural domain, with their delta method standard error:' );
ParamterNames = { 'rho_1'; 'rho_1Star'; 'psi_1'; 'psi_1Star'; 'sigma_1'; 'sigma_1Star'; 'sigma_infinity'; 'sigma_2' };
DeltaMethodStandardErrors = sqrt( diag( EstimatedParametersCovariance ) ) .* abs( DNaturalDomainParameters );
disp( table( NaturalDomainParameters, DeltaMethodStandardErrors, 'RowNames', ParamterNames ) );
disp( ' ' );

[ SmoothedState, ~, SmoothingOutput ] = smooth( EstimatedModel, Data, 'SquareRoot', true );

Pi = SmoothedState( :, 14 );
PiStar = SmoothedState( :, 15 );
VarPiStar = arrayfun( @( t ) SmoothingOutput( t ).SmoothedStatesCov( 15, 15 ), 1 : T, 'UniformOutput', true ).';
StdPiStar = sqrt( VarPiStar );

FillBetweenLines = @( X, Y1, Y2 ) fill( [ X(:).', flip( X(:) ).' ],  [ Y1(:).', flip( Y2(:) ).' ], [ 0.8, 0.8, 0.8 ], 'EdgeColor', 'none', 'LineStyle', 'none' );

figure;
FillBetweenLines( Dates, 12 * ( PiStar - norminv( 0.95 ) * StdPiStar ), 12 * ( PiStar + norminv( 0.95 ) * StdPiStar ) ); % 90% confidence bands.
hold on;
plot( Dates, 12 * Pi, 'k--', Dates, 12 * PiStar, 'k' );
xlabel( 'Year' );
ylabel( 'Annualized inflation rate' );

                % saveas( gcf, 'Outputs/SupplementalFigure7MonthlyPiPiStar.svg' );

%% Prepare to perform the main regression.

[ A, B, C, D, Mean0 ] = PiStarParameterFunction( EstimatedParameters );

FiveYearExpectation = SmoothedState( :, 16 ) + SmoothedState( :, 15 ); % Expectation at t of pi^*_{t-1} + pi^*_{t} +...

ForwardState = SmoothedState;

for t = 1 : ( 5 * 12 - 2 )
    ForwardState = ForwardState * A.';
    FiveYearExpectation = FiveYearExpectation + ForwardState( :, 15 );
end

FiveYearExpectation = FiveYearExpectation / 5;

assert( ismember( Dates( 1 ), CCRateT5YIEM.Properties.RowTimes ) );
BreakevenInflation = CCRateT5YIEM.CCRateT5YIEM( ismember( CCRateT5YIEM.Properties.RowTimes, Dates ) );

Data = readtable( 'Inputs/OilSupplyNewsShocksKaenzig.xlsx', 'Sheet', 'Monthly', 'ReadVariableNames', true, 'VariableNamingRule', 'modify' );

OilSupplyNewsShockDates = repmat( datetime( 1000, 1, 1 ), size( Data, 1 ), 1 );

for t = 1 : size( Data, 1 ) 
    
    Date = sscanf( Data.Date{ t }, '%dM%d' );
    assert( numel( Date ) == 2 );

    OilSupplyNewsShockDates( t ) = datetime( Date( 1 ), Date( 2 ), 1 );

end

assert( ismember( Dates( 1 ), OilSupplyNewsShockDates ) );
OilSupplyNewsShock = Data.OilSupplyNewsShock( ismember( OilSupplyNewsShockDates, Dates ) );

assert( ismember( Dates( 1 ), GS5.Properties.RowTimes ) );
TreasuryYield = GS5.CCRateGS5( ismember( GS5.Properties.RowTimes, Dates ) );

disp( 'December 2021 SEP inflation forecast for 2021:' );
disp( HistoricalSEPSRData( year( HistoricalSEPDates ) == 2021, year( HistoricalSEPVintages ) == 2021 & month( HistoricalSEPVintages ) == 12 ) );
disp( 'Realised 2021 inflation:' );
Indices2021 = find( year( Dates ) == 2021 );
disp( mean( sum( [ Pi( Indices2021 ), Pi( Indices2021 - 1 ), Pi( Indices2021 - 2 ) ] ) ) );
disp( 'March and June 2021 SEP inflation forecast for 2021:' );
disp( HistoricalSEPSRData( year( HistoricalSEPDates ) == 2021, year( HistoricalSEPVintages ) == 2021 & ( month( HistoricalSEPVintages ) == 3 | month( HistoricalSEPVintages ) == 6 ) ) );

% Now convert to quarterly.

AggregationVector = [ 1; 1; 1 ];

PiStar = SmoothedState( :, 15 : 17 ) * AggregationVector;

VarPiStar = arrayfun( @( t ) AggregationVector.' * SmoothingOutput( t ).SmoothedStatesCov( 15 : 17, 15 : 17 ) * AggregationVector, 1 : T, 'UniformOutput', true ).';

StdPiStar = sqrt( VarPiStar );

                % VarObsRelativeInflation = arrayfun( @( t ) [ 12 -12 ] * SmoothingOutput( t ).SmoothedStatesCov( [ 14, 15 ], [ 14, 15 ] ) * [ 12; -12 ], 1 : T, 'UniformOutput', true ).';
                % VarObsRelativeInflation = [ NaN; VarObsRelativeInflation( 1 : ( end - 1 ) ) ];
                % VarRelativeInflation = [ 12 -12 ] * EstimatedModel.Cov0( [ 14, 15 ], [ 14, 15 ] ) * [ 12; -12 ];
                % LambdaRelativeInflation = VarRelativeInflation ./ ( VarRelativeInflation + VarObsRelativeInflation );

QuarterlySums = @( z ) [ z, [ NaN; z( 1 : ( end - 1 ) ) ], [ NaN; NaN; z( 1 : ( end - 2 ) ) ] ] * AggregationVector;

Pi = QuarterlySums( Pi );
OilSupplyNewsShock = QuarterlySums( OilSupplyNewsShock );

Dates = Dates( 11 : 3 : end ) + calendarDuration( 0, 0, 14 ); % Middle of quarter date.

PiLong = Pi( 9 : 3 : end );
Pi = Pi( 12 : 3 : end ); % First element is Q4.
OilSupplyNewsShock = OilSupplyNewsShock( 12 : 3 : end );
DPiStar = diff( PiStar( 9 : 3 : end ) );
PiStarLong = PiStar( 9 : 3 : end );
PiStar = PiStar( 12 : 3 : end );
StdPiStar = StdPiStar( 12 : 3 : end );

FiveYearExpectation = FiveYearExpectation( 8 : 3 : end );
BreakevenInflation = BreakevenInflation( 9 : 3 : end );
TreasuryYield = TreasuryYield( 9 : 3 : end );

disp( 'Q1 and Q2 2021 CPI Inflation:' );
disp( 4 * Pi( year( Dates ) == 2021 & month( Dates ) <= 6 ).' );
disp( ' ' );

figure;
FillBetweenLines( Dates, 4 * ( PiStar - norminv( 0.95 ) * StdPiStar ), 4 * ( PiStar + norminv( 0.95 ) * StdPiStar ) ); % 90% confidence bands.
hold on;
plot( Dates, 4 * Pi, 'k--', Dates, 4 * PiStar, 'k' );
xlabel( 'Year' );
ylabel( 'Annualized inflation rate' );
saveas( gcf, 'Outputs/SupplementalFigure7PiPiStar.svg' );

%% Now perform the main regression.

disp( ' ' );
disp( '*** MAIN PAPER SECTION 6 & SUPPLEMENTAL APPENDIX J.5 ***' );
disp( ' ' );

figure;
plot( Dates( 1 : numel( BreakevenInflation( 2 : end ) ) ), BreakevenInflation( 2 : end ), 'k--', Dates( 1 : numel( FiveYearExpectation( 2 : end ) ) ), FiveYearExpectation( 2 : end ), 'k' );
xlabel( 'Year' );
ylabel( 'Average annual inflation rate' );
saveas( gcf, 'Outputs/SupplementalFigure8FiveYearExpectation.svg' );

                % save Outputs/Results AlphaBetaDates alpha beta Dates PiStar FiveYearExpectation;

                % writetimetable( timetable( Dates, PiStar, FiveYearExpectation( 2 : end ) ), 'Outputs/EstimatedPiStarAndFiveYearExpectation.xlsx', 'WriteMode', 'replacefile' );

assert( numel( Dates ) == numel( Pi ) );
assert( numel( Dates ) == numel( PiStar ) );
assert( numel( Dates ) + 1 == numel( PiLong ) );
assert( numel( Dates ) + 1 == numel( PiStarLong ) );
assert( numel( Dates ) + 1 == numel( FiveYearExpectation ) );

T = numel( BreakevenInflation ) - 1;

assert( numel( TreasuryYield ) - 1 == T );
assert( numel( OilSupplyNewsShock ) <= T );

if numel( Dates ) > T

    Dates = Dates( 1 : T );
    Pi = Pi( 1 : T );
    PiStar = PiStar( 1 : T );
    PiLong = PiLong( 1 : ( T + 1 ) );
    PiStarLong = PiStarLong( 1 : ( T + 1 ) );
    FiveYearExpectation = FiveYearExpectation( 1 : ( T + 1 ) );

end

RelativeFiveYearExpectation = BreakevenInflation - FiveYearExpectation;

x = Pi - PiStar;
xLong = PiLong - PiStarLong;
z = OilSupplyNewsShock;
y = 0.25 * ( RelativeFiveYearExpectation( 2 : end ) - RelativeFiveYearExpectation( 1 : ( end - 1 ) ) );
f = 0.25 * RelativeFiveYearExpectation - xLong / 20;
yMod = f( 2 : end ) - f( 1 : ( end - 1 ) );

yA = 0.25 * RelativeFiveYearExpectation( 2 : end );
yB = 0.25 * ( BreakevenInflation( 2 : end ) - BreakevenInflation( 1 : ( end - 1 ) ) );
yC = 0.25 * BreakevenInflation( 2 : end );
yD = 0.25 * ( TreasuryYield( 2 : end ) - TreasuryYield( 1 : ( end - 1 ) ) );
yE = 0.25 * TreasuryYield( 2 : end );

% Drop the first year to give the smoother time to "burn-in"
x( 1 : 4 ) = [];
y( 1 : 4 ) = [];
z( 1 : 4 ) = [];
yMod( 1 : 4 ) = [];
yA( 1 : 4 ) = [];
yB( 1 : 4 ) = [];
yC( 1 : 4 ) = [];
yD( 1 : 4 ) = [];
yE( 1 : 4 ) = [];
Dates( 1 : 4 ) = [];

% Drop the observations not covered by the news shock series.
T = numel( z );

x( ( T + 1 ) : end ) = [];
y( ( T + 1 ) : end ) = [];
yMod( ( T + 1 ) : end ) = [];
yA( ( T + 1 ) : end ) = [];
yB( ( T + 1 ) : end ) = [];
yC( ( T + 1 ) : end ) = [];
yD( ( T + 1 ) : end ) = [];
yE( ( T + 1 ) : end ) = [];
Dates( ( T + 1 ) : end ) = [];

disp( 'Sample for the main regression (not including the observation consumed by the lag) (Y, Q, Y, Q):' );
disp( [ year( Dates( 1 ) ), quarter( Dates( 1 ) ), year( Dates( end ) ), quarter( Dates( end ) ) ] );
disp( ' ' );

disp( 'Minimum annual yield on 5-year treasuries:' )
disp( min( DGS5.DGS5 ) );
disp( ' ' );

Years = unique( year( Dates ) );

Labels = repmat( Years(:).', 4, 1 );
Labels = cellstr( [ num2str( Labels(:) ) repmat( [ repmat( 'Q', 4, 1 ) num2str( ( 1 : 4 ).' ) ], numel( Years ), 1 ) ] );
Labels( 1 : 3 ) = [];
Labels( ( numel( y ) + 1 ) : end ) = [];

[ theta, ~, ~, Ptheta ] = hacIV( x, x, y );
[ thetaMod, ~, ~, PthetaMod ] = hacIV( x, x, yMod );

[ thetaIV, ~, ~, PthetaIV, ~, ~, ~, ~, F ] = hacIV( x, z, y );
[ ~, ~, ~, PthetaIV01 ] = hacIV( x, z, y - 0.1 * x );
[ thetaIVMod, ~, ~, PthetaIVMod, ~, ~, ~, ~, FMod ] = hacIV( x, z, yMod );

figure;
autocorrNoMean( yMod );
figure;
autocorrNoMean( x );

resid = y - x * theta;
residIV = y - x * thetaIV;
residMod = yMod - x * thetaMod;
residIVMod = yMod - x * thetaIVMod;

% Now the approach using a monetary shock instrument

MS = readmatrix( 'Inputs/FOMC_Bauer_Swanson.xlsx', 'Sheet', 'Monthly SVAR Data', 'Range', 'I431:I565' );
LastMS = numel( MS ) / 3;
MS = sum( reshape( MS, 3, LastMS ) ).';

figure; plot( 1 : LastMS, residIV( 1 : LastMS ), 'k', 1 : LastMS, MS, 'r' );
figure; plot( 1 : LastMS, residIVMod( 1 : LastMS ), 'k', 1 : LastMS, MS, 'r' );

FirstMS = 1;

disp( 'Sample for the Bauer-Swanson regression (not including the observation consumed by the lag) (Y, Q, Y, Q):' );
disp( [ year( Dates( FirstMS ) ), quarter( Dates( FirstMS ) ), year( Dates( LastMS ) ), quarter( Dates( LastMS ) ) ] );
disp( ' ' );

MS = MS( FirstMS : LastMS );
xMS = x( FirstMS : LastMS, : );
yMS = y( FirstMS : LastMS );
yMSMod = yMod( FirstMS : LastMS );

NMS = numel( MS );

XMS = [ xMS, yMS ];
XMSMod = [ xMS, yMSMod ];

[ ab, var_ab ] = hacIV( XMS, XMS, MS );

thetaMS = -ab( 1 ) / ab( 2 );

Grad = [ -1 / ab( 2 ); ab( 1 ) / ( ab( 2 ) * ab( 2 ) ) ];

var_thetaMS = Grad.' * var_ab * Grad;

t_thetaMS = thetaMS ./ sqrt( var_thetaMS );

PthetaMS = 2 * ( 1 - tcdf( abs( t_thetaMS ), NMS - 2 ) );


[ ab, var_ab ] = hacIV( XMSMod, XMSMod, MS );

thetaMSMod = -ab( 1 ) / ab( 2 );

Grad = [ -1 / ab( 2 ); ab( 1 ) / ( ab( 2 ) * ab( 2 ) ) ];

var_thetaMSMod = Grad.' * var_ab * Grad;

t_thetaMSMod = thetaMSMod ./ sqrt( var_thetaMSMod );

PthetaMSMod = 2 * ( 1 - tcdf( abs( t_thetaMSMod ), NMS - 2 ) );


residMS = y - x * thetaMS;
residMSMod = yMod - x * thetaMSMod;

RSS = sum( resid .^ 2 );
RSSIV = sum( residIV .^ 2 );
RSSMS = sum( residMS .^ 2 );
RSSMod = sum( residMod .^ 2 );
RSSIVMod = sum( residIVMod .^ 2 );
RSSMSMod = sum( residMSMod .^ 2 );

disp( 'Results from the unmodified quarterly exercise (this is not the main specification used in the paper, these results are only mentioned in footnote 41):' );
disp( 'Correlation of y with z, p-value below:' );
disp( corrNoMean( y, z ) );
disp( 'Estimated theta from the quarterly exercise (OLS, IV, MS):' );
disp( [ theta thetaIV thetaMS ] );
disp( 'HAC P-value from the quarterly exercise (OLS, IV, MS):' );
disp( [ Ptheta PthetaIV PthetaMS ] );
disp( 'HAC P-value of theta=0.1 from the quarterly exercise (IV):' );
disp( PthetaIV01 );
disp( 'First stage IV F value from the quarterly exercise:' );
disp( F );
disp( 'Correlation of residuals with the Bauer Swanson monetary shocks (OLS, IV, MS), p-values below:' );
disp( [ corrNoMean( resid( FirstMS : LastMS ), MS ), corrNoMean( residIV( FirstMS : LastMS ), MS ), corrNoMean( residMS( FirstMS : LastMS ), MS ) ] );
disp( 'Spearman correlation of residuals with the Bauer Swanson monetary shocks (OLS, IV, MS):' );
disp( [ corr( resid( FirstMS : LastMS ), MS, 'Type', 'Spearman' ), corr( residIV( FirstMS : LastMS ), MS, 'Type', 'Spearman' ), corr( residMS( FirstMS : LastMS ), MS, 'Type', 'Spearman' ) ] );
disp( 'R^2 of changes in relative breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSS ./ sum( y .* y ), 1 - RSSIV ./ sum( y .* y ), 1 - RSSMS ./ sum( y .* y ) ] );
disp( 'R^2 of relative breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSS ./ sum( yA .* yA ), 1 - RSSIV ./ sum( yA .* yA ), 1 - RSSMS ./ sum( yA .* yA ) ] );
disp( 'R^2 of changes in breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSS ./ sum( yB .* yB ), 1 - RSSIV ./ sum( yB .* yB ), 1 - RSSMS ./ sum( yB .* yB ) ] );
disp( 'R^2 of breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSS ./ sum( yC .* yC ), 1 - RSSIV ./ sum( yC .* yC ), 1 - RSSMS ./ sum( yC .* yC ) ] );
disp( 'R^2 of changes in yields from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSS ./ sum( yD .* yD ), 1 - RSSIV ./ sum( yD .* yD ), 1 - RSSMS ./ sum( yD .* yD ) ] );
disp( 'R^2 of yields from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSS ./ sum( yE .* yE ), 1 - RSSIV ./ sum( yE .* yE ), 1 - RSSMS ./ sum( yE .* yE ) ] );
disp( ' ' );

disp( 'Results from the modified quarterly exercise (this is the main specification used in the paper):' );
disp( 'Correlation of yMod with z, p-value below:' );
disp( corrNoMean( yMod, z ) );
disp( 'Estimated thetaMod from the quarterly exercise (OLS, IV, MS):' );
disp( [ thetaMod thetaIVMod thetaMSMod ] );
disp( 'HAC P-value from the quarterly exercise (OLS, IV, MS):' );
disp( [ PthetaMod PthetaIVMod PthetaMSMod ] );
disp( 'First stage IV F value from the quarterly exercise:' );
disp( FMod );
disp( 'Correlation of residuals with the Bauer Swanson monetary shocks (OLS, IV, MS), p-values below:' );
disp( [ corrNoMean( residMod( FirstMS : LastMS ), MS ), corrNoMean( residIVMod( FirstMS : LastMS ), MS ), corrNoMean( residMSMod( FirstMS : LastMS ), MS ) ] );
disp( 'Spearman correlation of residuals with the Bauer Swanson monetary shocks (OLS, IV, MS):' );
disp( [ corr( residMod( FirstMS : LastMS ), MS, 'Type', 'Spearman' ), corr( residIVMod( FirstMS : LastMS ), MS, 'Type', 'Spearman' ), corr( residMSMod( FirstMS : LastMS ), MS, 'Type', 'Spearman' ) ] );
disp( 'R^2 of changes in modified relative breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( yMod .* yMod ), 1 - RSSIVMod ./ sum( yMod .* yMod ), 1 - RSSMSMod ./ sum( yMod .* yMod ) ] );
disp( 'R^2 of changes in relative breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( y .* y ), 1 - RSSIVMod ./ sum( y .* y ), 1 - RSSMSMod ./ sum( y .* y ) ] );
disp( 'R^2 of relative breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( yA .* yA ), 1 - RSSIVMod ./ sum( yA .* yA ), 1 - RSSMSMod ./ sum( yA .* yA ) ] );
disp( 'R^2 of changes in breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( yB .* yB ), 1 - RSSIVMod ./ sum( yB .* yB ), 1 - RSSMSMod ./ sum( yB .* yB ) ] );
disp( 'R^2 of breakeven inflation from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( yC .* yC ), 1 - RSSIVMod ./ sum( yC .* yC ), 1 - RSSMSMod ./ sum( yC .* yC ) ] );
disp( 'R^2 of changes in yields from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( yD .* yD ), 1 - RSSIVMod ./ sum( yD .* yD ), 1 - RSSMSMod ./ sum( yD .* yD ) ] );
disp( 'R^2 of yields from the quarterly exercise (OLS, IV, MS):' );
disp( [ 1 - RSSMod ./ sum( yE .* yE ), 1 - RSSIVMod ./ sum( yE .* yE ), 1 - RSSMSMod ./ sum( yE .* yE ) ] );
disp( ' ' );

figure;
scatter( 4 * x, 4 * y, 12, 'filled', 'k' );
hold on;
XLim = xlim;
YLim = ylim;
plot( XLim, XLim * theta, 'k:', 'LineWidth', 1 );
plot( XLim, XLim * thetaIV, 'k-', 'LineWidth', 1 );
plot( XLim, XLim * thetaMS, 'k--', 'LineWidth', 1 );
OnLeft = ismember( Labels, { '2008Q4', '2011Q4', '2012Q4', '2021Q1', '2020Q2', '2022Q1' } );
OnRight = ismember( Labels, { '2009Q1', '2009Q3', '2009Q4', '2015Q2', '2020Q1', '2021Q2', '2022Q2' } );
text( 4 * x( OnRight ) + 4 * 0.05, 4 * y( OnRight ), Labels( OnRight ), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( 4 * x( OnLeft ) - 4 * 0.05, 4 * y( OnLeft ), Labels( OnLeft ), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle' );
xlim( XLim );
ylim( YLim );
                % xlabel( '$x_t=\pi_t-\pi_t^*$ \textsf{(continuously-compounded annual rate)}', 'Interpreter', 'latex' );
                % ylabel( '$y_t$ \textsf{(continuously-compounded annual rate)}', 'Interpreter', 'latex' );
set( gca, 'Box', 'on' );
                % saveas( gcf, 'Outputs/PaperFigure1Quarterly.svg' );

figure;
scatter( 4 * x, 4 * yMod, 12, 'filled', 'k' );
hold on;
XLim = xlim;
YLim = ylim;
plot( XLim, XLim * thetaMod, 'k:', 'LineWidth', 1 );
plot( XLim, XLim * thetaIVMod, 'k-', 'LineWidth', 1 );
plot( XLim, XLim * thetaMSMod, 'k--', 'LineWidth', 1 );
OnLeft = ismember( Labels, { '2008Q4', '2011Q4', '2012Q2', '2012Q4', '2020Q2', '2021Q1', '2021Q3', '2022Q1', '2022Q3', '2022Q4' } );
OnRight = ismember( Labels, { '2009Q1', '2009Q3', '2009Q4', '2011Q1', '2011Q2', '2012Q3', '2013Q2', '2015Q2', '2017Q3', '2020Q1', '2021Q2', '2022Q2' } );
text( 4 * x( OnRight ) + 4 * 0.05, 4 * yMod( OnRight ), Labels( OnRight ), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle' );
text( 4 * x( OnLeft ) - 4 * 0.05, 4 * yMod( OnLeft ), Labels( OnLeft ), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle' );
xlim( XLim );
ylim( YLim );
                % xlabel( '$x_t=\pi_t-\pi_t^*$ \textsf{(continuously-compounded annual rate)}', 'Interpreter', 'latex' );
                % ylabel( '$y_t$ \textsf{(continuously-compounded annual rate)}', 'Interpreter', 'latex' );
set( gca, 'Box', 'on' );
saveas( gcf, 'Outputs/PaperFigure1QuarterlyMod.svg' );

%% Clean up.

if ~isempty( DynarePath )
    addpath( [ DynarePath '/matlab' ] );
    delete private/dynareNoLog.m;
end

delete private/autocorrNoMean.m;

warning( 'on', 'MATLAB:DELETE:FileNotFound' );
warning( 'on', 'MATLAB:table:ModifiedAndSavedVarnames' );
warning( 'on', 'MATLAB:RMDIR:RemovedFromPath' );
warning( 'on', 'MATLAB:MKDIR:DirectoryExists' );

disp( ' ' );
disp( 'Run completed successfully. Elapsed time in seconds:' );
disp( toc );
disp( ' ' );

diary off;
