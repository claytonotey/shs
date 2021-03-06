#include "Params.h"
#include <iostream>

using namespace std;

ParamValue paramNULLInst;

size_t mapHasher(const std::unordered_map<std::string, ParamValue> &m)
{ 
  size_t ret = 0;
  for(std::unordered_map<std::string, ParamValue>::const_iterator i = m.begin(); i != m.end(); ++i) {
    ret ^= (std::hash<std::string>()(i->first) ^ (i->second.hash() << 1));
  }
  return ret;
}

size_t paramsHasher(const Params &m) 
{
  return mapHasher(m.params);
}

ParamValue Params :: parseParam(const std::string &s)
{
  try {
    Complex z = parseComplex(s);
    return ParamComplex(z);
  } catch(std::invalid_argument e) {
    try {
      int i = std::stoi(s);
      double f = std::stof(s);
      if(i == f) {
           return ParamInt(i);
      } else {
        return ParamFloat(f);
      }
    } catch(std::invalid_argument e) {
      return ParamString(s);
    }
  }
  return paramNULLInst;
}

Params :: Params() {}

ParamValue Params :: operator[](const std::string &s) const {
  ParamMap::type::const_iterator i = params.find(s);
  if(i != params.end()) {
    return (i->second);
  } else {
    return paramNULLInst;
  }
}

ParamValue &Params :: operator[](const std::string &s) {
  return params[s];
}

bool Params::operator==(const Params& other) const {
  return params == other.params;
}

bool Params :: operator!=(const Params& other) const {
  return !(*this == other);
}

Params :: Params(const string &s) 
{
  int pos0 = 0;
  int pos1;
  int pos2;
  do {
    pos1 = s.find_first_of("=", pos0);
    pos2 = s.find_first_of(",", pos1);
    if(pos1 == string::npos) {
      throw invalid_argument("exception parsing parameters");
    } else {
      string key = s.substr(pos0,pos1-pos0);
      string val = s.substr(pos1+1,pos2-pos1-1);
      (*this)[key] = Params::parseParam(val);
    }
    if(pos2 == string::npos) {
      break;
    }
    pos0 = pos2+1;
  } while(true);
}
