template < typename T >
void NumericalSolver( Mesh &mesh, SimulationInfo sim_info)
{
  T::solve( mesh, sim_info );
}
