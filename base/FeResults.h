#pragma once

class FeResults
{
public:
  virtual ~FeResults() = default;
  virtual void write() const = 0;
protected:
  FeResults() = default;
  FeResults(FeResults const &other) = delete;
  FeResults(FeResults &&other) = delete;
  FeResults &operator=(FeResults const &other) = delete;
};