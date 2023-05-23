function [ r, key ] = niederreiter2_generate ( dim_num, n, key, nbits )

%*****************************************************************************80
%
%% niederreiter2_generate() generates a set of Niederreiter values.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 March 2003
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer DIM_NUM, the spatial dimension.
%
%    integer N, the number of points desired.
%
%    integer key, a key for the random
%    number generator.
%
%  Output:
%
%    real R(DIM_NUM,N), the points.
%
%    integer key, a key for the random number generator.
%
  r = zeros(dim_num,n);

  for j = 1 : n
    [ r(1:dim_num,j), key ] = niederreiter2 ( dim_num, key, nbits );
  end

  return
end
