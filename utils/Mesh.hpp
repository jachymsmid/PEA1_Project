#pragma once

// should i use enum and switch-case?

#include <array>
#include <string>

template
<
  std::size_t rows,
  std::size_t cols,
  class RealNumber,
  class IndexingType,
  class SimulationInfo
>
class Mesh
{
private:

  // would a array of 3 big arrays be better? or are 3 separate arrays better?
  std::array< RealNumber, rows * cols > data;
  std::array< RealNumber, rows > x_cord;
  std::array< RealNumber, cols > y_cord;

public:

  Mesh() = default;

  // find a better way?
  Mesh( SimulationInfo sim_info ){} 

  // is the copied mesh gonna have the same templates?
  Mesh( const Mesh &mesh ) : data(mesh.data),
          x_cord(mesh.x_cord),
          y_cord(mesh.y_cord) {}

  RealNumber& value_ref( size_t i, size_t j )
  {
    data[ IndexingType::index( i, j ) ];
  }

  const RealNumber &value( size_t i, size_t j ) const
  {
    data[ IndexingType::index( i, j ) ];
  }

  // coordinates are seemingly 2D but only the first, or the second index is used
  const RealNumber &x_coordinate( size_t i, size_t j ) const
  {
    x_cord[ IndexingType::index( i, j ) ];
  }

  const RealNumber &y_coordinate( size_t i, size_t j ) const
  {
    y_cord[ IndexingType::index( i, j ) ];
  }

  size_t getCols() const { return cols; }
  size_t getRows() const { return rows; }

  void construct_regular_grid( RealNumber x_step, RealNumber y_step )
  {
    for ( size_t i = 0; i < rows; i++ )
    {
      x_cord( IndexingType::index( i, 0 ) ) = i * x_step;
    }
    for ( size_t i = 0; i < cols; i++ )
    {
      y_cord( IndexingType::index( 0, i ) ) = i * y_step;
    }
  }

  void write_data( std::string file_name );
  // rewrite without using the VTK library

};



struct Column_major
{
public:
  static size_t index( size_t i, size_t j )
  {
    return Mesh::getCols() * i + j ;
  }
};

struct Row_major
{
};

struct Z_order 
{
};

