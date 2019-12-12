#pragma once

class SimulationResults
{
public:
  virtual ~SimulationResults() = default;
  virtual void write() const = 0;
protected:
  SimulationResults() = default;
  SimulationResults(SimulationResults const &other) = delete;
  SimulationResults(SimulationResults &&other) = delete;
  SimulationResults &operator=(SimulationResults const &other) = delete;
};