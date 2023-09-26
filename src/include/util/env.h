#pragma once

#include <string>
#include <optional>
#include <stdexcept>

namespace util {
    std::string getenv_exc(std::string name);

    std::string getenv_exc(std::string name) {
      char* value = std::getenv(name.c_str());

      if (value == nullptr) {
        throw std::runtime_error("Environment variable " + name + " not found!");
      } else {
        return std::string(value);
      }
    }
}