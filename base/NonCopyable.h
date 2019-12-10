#pragma once

class NonCopyable
{
public:
  virtual ~NonCopyable() = default;

protected:
  NonCopyable(NonCopyable const &) = delete;
  NonCopyable &operator=(NonCopyable const &) = delete;

  NonCopyable(NonCopyable &&) = delete;
  NonCopyable &operator=(NonCopyable &&) = delete;
};