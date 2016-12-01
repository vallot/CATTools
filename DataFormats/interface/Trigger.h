#ifndef CATTools_TriggerBits_H
#define CATTools_TriggerBits_H

#include <vector>
#include <string>
#include <map>

namespace cat {

typedef std::map<std::string, unsigned short> TriggerResValues;

class TriggerNames
{
public:
  TriggerNames();
  TriggerNames(const TriggerResValues& results) { set(results); }
  virtual ~TriggerNames() {};

  // Getters
  int index(const std::string& name) const;
  std::string name(const size_t i) const { return names_.at(i); }
  size_t size() const { return names_.size(); }
  std::vector<std::string> names() const { return names_; }

  int print() const;

  // Setters
  void set(const TriggerResValues& results);

private:
  std::vector<std::string> names_;
};

class TriggerBits
{
public:
  TriggerBits();
  TriggerBits(const TriggerResValues& results) { set(results); }
  virtual ~TriggerBits() {};

  // Setters
  void set(const TriggerResValues& results);

  // Getters
  unsigned short result(const int i) const { return values_.at(i); } 

private:
  std::vector<unsigned short> values_;

};

}

#endif
