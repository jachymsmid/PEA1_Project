template < typename T >
void BoundaryConditions( Mesh &mesh )
{
  T::impose( mesh );
};
