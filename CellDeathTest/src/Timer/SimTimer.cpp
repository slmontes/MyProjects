#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include "SimulationTime.hpp"
#include "SimTimer.hpp"

void SimTimer::MoveOneTimeStep()
{
  SimulationTime* p_simulation_time = SimulationTime :: Instance();
  p_simulation_time->IncrementTimeOneStep();
  return;
}
int SimTimer::IsSimFinished(int sim_ongoing)
{
  SimulationTime* p_simulation_time = SimulationTime :: Instance();
  if (p_simulation_time->IsFinished())
  {
    sim_ongoing=0;
  }
  return sim_ongoing;
}
