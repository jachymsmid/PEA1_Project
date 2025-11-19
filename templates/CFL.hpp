template < typename T >
RealNumber CFL( RealNumber dx, RealNumber dy )
{
  return T::cfl( dx, dy );
}
