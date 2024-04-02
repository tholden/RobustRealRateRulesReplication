function [ rho_p ] = corrNoMean( x, y )

    n = size( x, 1 );
    assert( n == size( y, 1 ) );
    
    assert( size( x, 2 ) == size( y, 2 ) );
    
    rho = sum( x .* y ) ./ sqrt( sum( x .* x ) .* sum( y .* y ) );
    
    t = rho .* sqrt( ( n - 2 ) ./ ( 1 - rho .^ 2 ) );
    p = 2 * tcdf( -abs( t ), n - 2 );
    
    rho_p = [ rho; p ];

end
