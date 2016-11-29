#include <cassert>

#include "profile_h/Registers.h"

Registers::~Registers()
{
  for (auto it = regs.begin(); it != regs.end(); ++it)
  {
    delete it->second;
  }
}

void Registers::createRegister(
    std::string baseName, unsigned int num_bytes)
{
  assert(regs.find(baseName) == regs.end());
  unsigned int new_id = regs.size();
  regs[baseName] = new Register(baseName, num_bytes, new_id);
}

void Registers::getRegisterNames(std::vector<std::string> &names)
{
  for (auto it = regs.begin(); it != regs.end(); it++)
    names.push_back(it->first);
}

bool Registers::has(std::string baseName)
{
  return (regs.find(baseName) != regs.end());
}
