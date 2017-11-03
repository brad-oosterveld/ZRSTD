#!/bin/python

import math

#REV: Last K rows take up the last K*K/2? I.e. number of elements is
#for size K, not K*K/2 (which would be only half of middle row, but rather,
#K*(K-1)/2? Yes.

#REV: OK, given index i, there are K(K-1)/2-1-i guys greater than me. This is
#i' We want to find largest k, s.t. k(k-1)/2 <= i'

#http://stackoverflow.com/questions/242711/algorithm-for-index-numbers-of-triangular-matrix-coefficients

#k(k-1) <= 2*ii means (k*k - k) <= 2*ii, err, (k-(1/2))^2 - (1/4)?
#Because, um, k-(1/2) * k-(1/2) = k^2 - k/2 - k/2 + 1/4.
#Which is, k^2 - k + 1/4. Yea so need to subtract...? -k diff I think.
#So, (k-(1/2))^2 - (1/4) <= 2*ii, or...
#k-(1/2) <= sqrt( 2*ii + 1/4 )
#k <= sqrt((2*ii + 1/4)*4 * (1/4)) + (1/2)
#k <= sqrt((8*ii + 1))/2 + 1/2
#k <= (sqrt(8*ii + 1) + 1)/2

#So, row index(i, M) is M-1-largest(k)


#for column, K(K-1)/2 elements out of ii are in later rows, so
# jj=ii-K(K-1)/2 later elements in "this" row... (which has K elements..)
# so, K-jj;


def rowidx( i, M ):
    m = float(M);
    ii = ( m*( m-1.0 ))/2.0 - 1.0 - float(i);
    K = math.floor( (math.sqrt(8.0*ii + 1) + 1)/2.0 );
    return int(M-1-K);
    
#    m = float(M);
#    row = ( (-2.0*m) - 1 + math.sqrt( ( (4.0*m*(m+1)) - (8.0*float(i)) - 7.0) ) ) / -2.0;
#    if( row == float( int( row ) ) ):
#        row = row - 1;
#    return int(row);

def colidx( i, M ):
    m = float(M);
    ii = float(m*(m-1)/2.0 - 1.0 - float(i));
    K = math.floor( (math.sqrt(8.0*ii + 1) + 1)/2.0);
    jj = int(ii) - (K*(K-1))/2.0;
    return int(K-jj-1);
#    row = rowidx( i, M );
#    col = i - ( M * row ) + (row*(row)) / 2;
#    return col;


#REV: Test a few:

N=3;
print( "For N=3: idx 2, (should be 1, 0): row=", rowidx( 2, N ), " col=", colidx( 2, N ) );

print( "For N=3: idx 1, (should be 0, 1): row=", rowidx( 1, N ), " col=", colidx( 1, N ) );


N=4;
print( "For N=4: idx 2, (should be 0, 2): row=", rowidx( 2, N ), " col=", colidx( 2, N ) );

print( "For N=4: idx 4, (should be 1, 1): row=", rowidx( 4, N ), " col=", colidx( 4, N ) );

print( "For N=4: idx 5, (should be 2, 0): row=", rowidx( 5, N ), " col=", colidx( 5, N ) );

