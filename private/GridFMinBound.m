function [ x, fx ] = GridFMinBound( f, LB, UB, n, Verbosity, TolX )

    xs  = linspace( LB, UB, n );
    fxs = Inf( 1, n );
    
    if Verbosity > 0
        
        disp( ' ' );
        disp( 'Initial grid search:' );
        disp( ' ' );
        
    end
    
    for i = 2 : ( n - 1 )
        
        fxs( i ) = f( xs( i ) );
        
        if Verbosity > 0
            disp( [ xs( i ) fxs( i ) ] );
        end
        
    end
    
    fx = min( fxs );
    LB = xs( find( fxs == fx, 1, 'first' ) - 1 );
    UB = xs( find( fxs == fx, 1, 'last'  ) + 1 );
    
    if Verbosity > 0
        fminbndOptions = optimset( 'Display', 'iter', 'TolX', TolX );
    else
        fminbndOptions = optimset( 'Display', 'off', 'TolX', TolX );
    end
    
    if Verbosity > 0
        
        disp( ' ' );
        disp( 'Fine search:' );
        disp( [ LB UB ] );
        
    end
    
    [ x, fx ] = fminbnd( f, LB, UB, fminbndOptions );
    
end
