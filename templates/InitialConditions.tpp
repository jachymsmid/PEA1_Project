template < typename T >
void InitialConditions( Mesh &mesh )
{
  T::impose( mesh );
}
