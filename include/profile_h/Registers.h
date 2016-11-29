#ifndef __REGISTERS_H__
#define __REGISTERS_H__

#include <map>
#include <string>
#include <vector>

class Register
{
  public:
    Register(std::string baseName, unsigned int size, unsigned int id)
    {
      this->baseName = baseName;
      this->id = id;
      num_bytes = size;
      loads = 0;
      stores = 0;
    }

    void increment_loads() { loads++; }
    void increment_stores() { stores++; }

  private:
    std::string baseName;
    unsigned int num_bytes;
    unsigned int id;

    unsigned int loads;
    unsigned int stores;
};

class Registers
{

  public:
    Registers() {}
    ~Registers();

    void createRegister(
        std::string baseName, unsigned num_bytes);
    void getRegisterNames(std::vector<std::string> &names);
    Register* getRegister(std::string baseName) { return regs[baseName]; }
    bool has(std::string baseName);

  private:
    std::map<std::string, Register*> regs;

};

#endif
