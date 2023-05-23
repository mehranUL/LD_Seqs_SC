Reference: https://people.sc.fsu.edu/~jburkardt/m_src/niederreiter2/niederreiter2.html

function [ quasi, key_new ] = niederreiter2 ( dim, key, nbits )

%*****************************************************************************80
%
%% niederreiter2() returns an element of the Niederreiter sequence base 2.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    27 February 2014
%
%  Author:
%
%    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
%    MATLAB version by John Burkardt;
%    Performance enhancements by Jeremy Dewar, Tulane University.
%
%  Reference:
%
%    Harald Niederreiter,
%    Low-discrepancy and low-dispersion sequences,
%    Journal of Number Theory,
%    Volume 30, 1988, pages 51-70.
%
%  Input:
%
%    integer DIM, the dimension of the sequence to be generated.
%
%    integer key, the index of the element to compute.
%
%  Output:
%
%    real QUASI(DIM), the next quasirandom vector.
%
%    integer key_NEW, the next value of the key.
%
%  Local:
%
%    integer MAXDIM, the maximum dimension that will be used.
%
%    integer NBITS, the number of bits (not counting the sign) in a
%    fixed-point integer.
%
%    real RECIP is the multiplier which changes the
%    integers in NEXTQ into the required real values in QUASI.
%
%    integer NR_cj(MAXDIM,NBITS), the packed values of
%    Niederreiter's C(I,J,R).
%
%    integer NR_dim, the spatial dimension of the sequence
%    as specified on an initialization call.
%
%    integer NR_nextq(MAXDIM), the numerators of the next item in the
%    series.  These are like Niederreiter's XI(N) (page 54) except that
%    N is implicit, and the NR_nextq are integers.  To obtain
%    the values of XI(N), multiply by RECIP.
%
  global NR_cj;
  global NR_dim;
  global NR_nextq;
  global NR_key;

  maxdim = 20;
  %nbits = 27; % default value:31
  recip = 2.0^(-nbits);
%
%  Initialization.
%
  if ( isempty ( NR_dim ) || dim ~= NR_dim || key <= 0 )

    if ( dim <= 0 || maxdim < dim )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'NIEDERREITER2 - Fatal error!\n' );
      fprintf ( 1, '  Bad spatial dimension.\n' );
      error ( 'Fatal error in NIEDERREITER2!' );
    end

    NR_dim = dim;

    if ( key < 0 )
      key = 0;
    end

    NR_key = key;
%
%  Calculate the C array.
%
    NR_cj(1:dim,1:nbits) = calcc2 ( dim );

  end
%
%  Set up NEXTQ appropriately, depending on the Gray code of key.
%
%  You can do this every time, starting NEXTQ back at 0,
%  or you can do it once, and then carry the value of NEXTQ
%  around from the previous computation.
%
  if ( key ~= NR_key + 1 )

    gray = bitxor ( floor ( key ), floor ( key / 2 ) );

    NR_nextq(1:dim) = 0;

    r = 0;

    while ( gray ~= 0 )

      if ( rem ( gray, 2 ) ~= 0 )
        for i = 1 : dim
          NR_nextq(i) = bitxor ( floor ( NR_nextq(i) ), floor ( NR_cj(i,r+1) ) );
        end
      end

      gray = floor ( gray / 2 );
      r = r + 1;

    end

  end
%
%  Multiply the numerators in NEXTQ by RECIP to get the next
%  quasi-random vector.
%
  quasi(1:dim) = NR_nextq(1:dim) * recip;
%
%  Find the position of the right-hand zero in key.  This
%  is the bit that changes in the Gray-code representation as
%  we go from key to key + 1.
%
  r = 0;
  i = key;

  while ( rem ( i, 2 ) ~= 0 )
    r = r + 1;
    i = floor ( i / 2 );
  end
%
%  Check that we have not passed 2**NBITS calls.
%
  if ( nbits <= r )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'NIEDERREITER2 - Fatal error!\n' );
    fprintf ( 1, '  Too many calls!\n' );
    error ( 'Fatal error in NIEDERREITER2!' );
  end
%
%  Compute the new numerators in vector NEXTQ.
%
  for i = 1 : dim
    NR_nextq(i) = bitxor ( floor ( NR_nextq(i) ), floor ( NR_cj(i,r+1) ) );
  end

  NR_key = key;
  key_new = key + 1;

  return
end
