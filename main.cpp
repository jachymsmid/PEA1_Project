#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <sstream>

// TODO:  - [ ] better mesh representation
//            - [x] change the array of structer paradigm to structure of array
//                - [ ] how to acces the data?
//            - [x] copy constructor
//            - [x] implement flattened arrays 
//               - [x] maybe Z-curves?
//                  - [ ] conditions for Z-curves usage
//            - [ ] irregular grid method - too elaborate?
//        - [ ] upwind method
//        - [x] lax-friedrichs method
//        - [ ] lax-wendroff method
//            - [ ] calculate the spatial step in the function - preparation for irregular grid
//            - [ ] clearer and simpler code
//        - [x] SimulationInfo struct
//            - [x] constructor - is it working??
//        - [x] how to implement the numerical solver
//            - [x] using the function template
//            - [o] using a pointer to a function
//        - [x] how to implement boundary/initial conditions?
//            - [x] again function templates
//        - [ ] VTK library for data storage and visualization
//        - [ ] error handling
//        - [ ] generalize the schemes?

using RealNumber = float;

// choose discretization scheme - options are:
//      - Lax_Wendroff
//      - Lax_Friedrich
//      - Upwind

// struct to hold infrmation about the simulation
struct SimulationInfo
{
  RealNumber step_t, step_x, step_y, x_min, y_min, x_max, y_max;
  int number_x, number_y;

  // constructor
  SimulationInfo( const RealNumber step_t,
                  const RealNumber step_x,
                  const RealNumber step_y,
                  const RealNumber x_min,
                  const RealNumber y_min,
                  const RealNumber x_max,
                  const RealNumber y_max,
                  const int number_x,
                  const int number_y)
                : step_t(step_t),
                  step_x(step_x),
                  step_y(step_y),
                  x_min(x_min),
                  y_min(y_min),
                  x_max(x_max),
                  y_max(y_max),
                  number_x(number_x),
                  number_y(number_y)
  {
    std::cout << "Simulation info initialization..." << std::endl;
  }
};

// ---------------
//    mesh class
// ---------------
class Mesh
{
private:

    std::vector< RealNumber > data;
    std::vector< RealNumber > x_cord;
    std::vector< RealNumber > y_cord;
    size_t rows, cols;

    size_t indices(size_t n)
    {
      size_t result = 0;
      size_t base = 1;
      while (n > 0) 
      {
        if (n & 1) result += base;
        base *= 4;
        n >>= 1;
      }
      return result;
    }

public:
    // constructor
    Mesh(SimulationInfo sim_info)
        : data(sim_info.number_x * sim_info.number_y),
          x_cord(sim_info.number_x * sim_info.number_y),
          y_cord(sim_info.number_x * sim_info.number_y),
          rows(sim_info.number_x), cols(sim_info.number_y)
  {
    std::cout << "Mesh initialization..." << std::endl;
  }

    // copy constructor
    Mesh(Mesh &mesh)
        : data(mesh.data),
          x_cord(mesh.x_cord),
          y_cord(mesh.y_cord),
          rows(mesh.rows),
          cols(mesh.cols) {}

    // getter/setter using Z-order curve
    //RealNumber& value_ref(size_t i, size_t j)
    //{
    //  if ( 2 * indices(i) + indices(j) > rows * cols )
    //  {
    //    throw std::out_of_range("Index out of bounds");
    //  }
    //  return data[ 2 * indices(i) + indices(j) ];
    //}

    //const RealNumber& value( size_t i, size_t j )
    //{
    //if ( 2 * indices(i) + indices(j) > rows * cols )
    //  {
    //    throw std::out_of_range("Index out of bounds");
    //  }
    //  return data[ 2 * indices(i) + indices(j) ];
    //}

    //const RealNumber& x(size_t i, size_t j)
    //{
    //  if ( 2 * indices(i) + indices(j) > rows * cols )
    //  {
    //    throw std::out_of_range("Index out of bounds");
    //  }
    //  return x_cord[ 2 * indices(i) + indices(j) ];
    //}

    //const RealNumber& y(size_t i, size_t j)
    //{
    //  if ( 2 * indices(i) + indices(j) > rows * cols )
    //  {
    //    throw std::out_of_range("Index out of bounds");
    //  }
    //  return data[ 2 * indices(i) + indices(j) ];
    //}

  // getter/setter for data
    RealNumber& value_ref(size_t i, size_t j)
    {
      if ( i > rows || j > cols )
      {
        throw std::out_of_range("Index out of bounds");
      }

      return data[i * cols + j];
    }

    const RealNumber& value(size_t i, size_t j) const
    {
      if ( i > rows || j > cols )
      {
        throw std::out_of_range("Index out of bounds");
      }

      return data[i * cols + j];
    }

    // getter for x coordinates
    const RealNumber& x(size_t i, size_t j) const
    {
      if ( i >= rows || j >= cols )
      {
        throw std::out_of_range("Index out of bounds");
      }

      return x_cord[i * cols + j];
    }

    // getter for y coordinates
    const RealNumber& y(size_t i, size_t j) const
    {
      if ( i >= rows || j >= cols )
      {
        throw std::out_of_range("Index out of bounds");
      }
      return data[i * cols + j];
    }

    // getters for number of columns and rows
    size_t getCols() const { return cols; }
    size_t getRows() const { return rows; }

    // method for regular grid construction
    void construct_regular_grid(SimulationInfo sim_info)
    {
      for (size_t i = 0; i < rows; ++i)
      {
        for (size_t j = 0; j < cols; ++j)
        {
            x_cord[i * cols + j] = sim_info.x_min + i * sim_info.step_x;
            y_cord[i * cols + j] = sim_info.y_min + j * sim_info.step_y;
        }
      }
      std::cout << "Regular grid constructed..." << std::endl;
    }

    // write the data to a file
    void write_data(std::string file_name)
    {
      std::ofstream file;
      file.open(file_name);

      if (file)
      {
        for (size_t i = 0; i < rows; i++)
        {
          for (size_t j = 0; j < cols; j++)
          {
            file << x_cord[i * cols + j] << "," << y_cord[i * cols + j] << "," << data[i * cols + j] << std::endl;
          }
        }
      }
      file.close();
    }

    // print the mesh to the command line, for troubleshooting
    void print_mesh()
    {
      for (size_t i = 0; i < rows; i++)
      {
        for (size_t j = 0; j < cols; j++)
        {
          std::cout << data[ i * cols + j ] << " ";
        }
        std::cout << std::endl;
      }
    }
};

// ----------------------------
// numerical solvers structures
// ----------------------------

// lax-friedrichs scheme
struct Lax_Friedrichs
{
    static void solve(Mesh &mesh, SimulationInfo sim_info)
    {
        Mesh U_n(mesh); 
        RealNumber dx, dy, dt, p = 1.0f;

        dx = sim_info.step_x;
        dy = sim_info.step_y;
        dt = sim_info.step_t;

        for(size_t i = 1; i < mesh.getRows() - 1; i++)
        {
            for(size_t j = 1; j < mesh.getCols() - 1; j++)
            {

                RealNumber x_j = mesh.x(i, j); 

                mesh.value_ref(i,j) = 
                      1.0f / 4.0f * ( U_n.value(i + 1, j) + U_n.value(i - 1, j) + U_n.value(i, j + 1) + U_n.value(i, j - 1))
                    - dt * (( U_n.value(i + 1, j) - U_n.value(i - 1, j)) / ( 2.0f * dx ) // u_x
                        - 2.0f * p * x_j * ( U_n.value(i, j + 1) - U_n.value(i, j - 1)) / ( 2.0f * dy ) 
                      );
            }
        }
    }

    // Funkce pro výpočet maximálního stabilního časového kroku (dt)
    static RealNumber cfl(RealNumber dx, RealNumber dy)
    {
        // podmínka pozitivity
        RealNumber dt_limit_x = dx / 2.0f;
        RealNumber dt_limit_y = dy / 6.0f; 
        
        return std::min(dt_limit_x, dt_limit_y);
    }
};

// lax-wendroff scheme
struct Lax_Wendroff
{
  static void solve(Mesh &mesh, SimulationInfo sim_info)
  {
    Mesh U_n(mesh);
    RealNumber lx, ly, dx, dy, dt;
    dx = sim_info.step_x;
    dy = sim_info.step_y;
    dt = sim_info.step_t;
    lx = dt / dx;

    for(size_t i = 1; i < mesh.getRows() - 1; i++)
    {
      for(size_t j = 1; j < mesh.getCols() - 1; j++)
      {
        ly = 2*mesh.x(i,j)*dt/dy;
        mesh.value_ref(i,j) = U_n.value(i,j)
                - lx / 2 * ( U_n.value(i + 1,j) - U_n.value(i - 1,j) )
                - ly / 2 * ( U_n.value(i,j + 1) - U_n.value(i,j - 1) )
                + pow(lx,2.0) / 2 * ( U_n.value(i + 1,j) - 2 * U_n.value(i,j) + U_n.value(i - 1,j) )
                + pow(ly,2.0) / 2 * (U_n.value(i,j + 1) - 2 * U_n.value(i,j) + U_n.value(i,j - 1) )
                + ly * lx / 4 * (U_n.value(i + 1,j + 1) - U_n.value(i - 1,j + 1) - U_n.value(i + 1,j - 1) + U_n.value(i - 1,j - 1) );
      }
    }
  }

  static RealNumber cfl(RealNumber dx, RealNumber dy)
  {
    return 1.0 / (1.0 / dx + 1.0 / 2.0 * 1.5 * dy);
  }
};

// upwind scheme
struct Upwind
{
    static void solve(Mesh &mesh, SimulationInfo sim_info)
    {
        Mesh U_n(mesh);  // uložení stavu z předchozího kroku
        RealNumber dx = sim_info.step_x;
        RealNumber dy = sim_info.step_y;
        RealNumber dt = sim_info.step_t;

        // rychlosti šíření
        RealNumber a = 1.0f;   // směr proudění ve směru x
        RealNumber b = 1.0f;   // směr proudění ve směru y

        for (size_t i = 1; i < mesh.getRows() - 1; i++)
        {
            for (size_t j = 1; j < mesh.getCols() - 1; j++)
            {
                mesh.value_ref(i, j) =
                    U_n.value(i, j)
                    - a * (dt / dx) * (U_n.value(i, j) - U_n.value(i - 1, j))
                    - b * (dt / dy) * (U_n.value(i, j) - U_n.value(i, j - 1));
            }
        }
    }

    static RealNumber cfl(RealNumber dx, RealNumber dy)
    {
        RealNumber a = 1.0f, b = 1.0f;
        //podmínka: dt <= 1 / ( |a|/dx + |b|/dy )
        return 1.0f / (fabs(a)/dx + fabs(b)/dy);
    }
};

// general function to call specific numerical scheme
template <typename T>
void NumericalSolver(Mesh &mesh, SimulationInfo sim_info)
{
  T::solve(mesh, sim_info);
}

using NumericalScheme = Lax_Wendroff;


// generalized cfl condition function
template <typename T>
RealNumber CFL( RealNumber dx, RealNumber dy)
{
  return T::cfl(dx, dy);
}

// -----------------------------
// initial conditions structures
// -----------------------------
struct My_Initial_Conditions
{
  static void impose(Mesh &mesh)
  {
    for (size_t i = 0; i < mesh.getRows(); i++)
    {
      for (size_t j = 0; j < mesh.getCols(); j++)
      {
        mesh.value_ref(i,j) = 100 * exp( - ( pow(mesh.x(i,j) + 1, 2.0) + pow(mesh.y(i,j), 2.0) ) / 0.01);
      }
    }
  }
};

// generalized initial conditions function
template< typename T >
void InitialConditions( Mesh &mesh)
{
  T::impose( mesh );
}

// ------------------------------
// boundary conditions structures
// ------------------------------
struct Zeros 
{
  static void impose(Mesh &mesh)
  {
    for (size_t i = 0; i < mesh.getRows(); ++i)
    {
      for (size_t j = 0; j < mesh.getCols(); ++j)
      {
        if (i == 0 || i == mesh.getRows() - 1 || j == 0 || j == mesh.getCols() - 1)
        {
          mesh.value_ref(i,j) = 0; 
        }
      }
    }
  }
};

// generalized boundary conditions function 
template< typename T >
void BoundaryConditions( Mesh &mesh )
{
  T::impose( mesh );
}

// -------------
//      main
// -------------
int main()
{
  // domain definition
  const RealNumber x_min = - 1.5, x_max = 1.5, y_min = - 1.5, y_max = 1.5;
  // number of grid points in the x, y direction respectively
  const int N_x = 64, N_y = 64;
  const RealNumber dx = (x_max-x_min) / N_x, dy = (y_max-y_min) / N_y;
  const RealNumber T = 2.f;
  RealNumber t = 0.f;
  const std::string sim_name;
  RealNumber dt = CFL < NumericalScheme >(dx, dy);

  SimulationInfo sim_info(dt, dx, dy, x_min, y_min, x_max, y_max, N_x, N_y);
  

  // initialize the mesh, U and file_name
  Mesh mesh(sim_info); 
  std::string file_name;

  mesh.construct_regular_grid(sim_info);

  InitialConditions< My_Initial_Conditions >(mesh);
  std::cout << "Imposing initial conditions..." << std::endl;

  mesh.write_data("sim/output_t_0.00000.csv");

  // std::cout << "Enter simulation name (folder with the same name will be created): ";
  // std::cin << sim_name;

  // main loop
  while (t <= T)
  {
    t += dt;
    NumericalSolver< NumericalScheme >(mesh, sim_info);
    BoundaryConditions< Zeros >(mesh);

    // this is kinda awkward
    std::ostringstream fn;
    fn << "sim/output_t_" << std::fixed << std::setprecision(5) << t << ".csv";
    file_name = fn.str();

    mesh.write_data(file_name);
  }
  return 0;
}
