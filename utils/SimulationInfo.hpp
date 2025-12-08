#pragma once
template < class RealNumber >
struct SimulationInfo
{
  const RealNumber step_t;
  const RealNumber step_x;
  const RealNumber step_y;
  const RealNumber x_min;
  const RealNumber y_min;
  const RealNumber x_max;
  const RealNumber y_max;
  const int number_x;
  const int number_y;
  const RealNumber x_speed;
  const RealNumber y_speed;

  SimulationInfo( RealNumber step_t,
                  RealNumber step_x,
                  RealNumber step_y,
                  RealNumber x_min,
                  RealNumber y_min,
                  RealNumber x_max,
                  RealNumber y_max,
                  int number_x,
                  int number_y,
                  RealNumber x_speed,
                  RealNumber y_speed)
                : step_t(step_t),
                  step_x(step_x),
                  step_y(step_y),
                  x_min(x_min),
                  y_min(y_min),
                  x_max(x_max),
                  y_max(y_max),
                  number_x(number_x),
                  number_y(number_y),
                  x_speed(x_speed),
                  y_speed(y_speed) {}
};
