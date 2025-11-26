#include "GeneralVector.hpp"
#include <string>
#include <vector>
#include "SimulationInfo.hpp"


template < class RealNumber >
class Mesh
{
private:

  std::vector< RealNumber > data;
  std::vector< RealNumber > x_cord;
  std::vector< RealNumber > y_cord;
  std::size_t rows, cols;

public:

  Mesh( SimulationInfo< RealNumber > sim_info )
    : data(sim_info.number_x * sim_info.number_y),
      x_cord(sim_info.number_x * sim_info.number_y),
      y_cord(sim_info.number_x * sim_info.number_y),
      rows(sim_info.number_x),
      cols(sim_info.number_y) {}

  Mesh( const Mesh &mesh )
    : data(mesh.data),
      x_cord(mesh.x_cord),
      y_cord(mesh.y_cord),
      rows(mesh.rows),
      cols(mesh.cols) {}

  RealNumber &value_ref( std::size_t i, std::size_t j ) { return data[i * cols + j]; }

  const RealNumber &value( size_t i, size_t j ) const { return data[i*cols + j]; }

  const RealNumber &x( size_t i, size_t j ) const { return x_cord[i*cols + j]; }

  const RealNumber &y( size_t i, size_t j ) const { return y_cord[i*cols + j]; }

  size_t getCols() const { return cols; }
  size_t getRows() const { return rows; }

  void construct_regular_grid( SimulationInfo< RealNumber > sim_info )
  {
    for ( size_t i = 0; i < rows * cols; i++)
    {
    }
  }

  void write_data( std::string file_name );

};
