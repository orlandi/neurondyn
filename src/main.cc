/*
 * Copyright (c) 2009-2015 Javier G. Orlandi <javierorlandi@javierorlandi.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include "main.h"
#include "netdyn.h"

#include <cassert>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdio.h>
 
int main(int argc, char* argv[])
{
  static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");
  // Disable buffering
  setvbuf(stdout, NULL, _IOLBF, 0);
  int i = 1;
  std::stringstream seedsUsed, configFile;
  while (i > 0)
  {
    simulation = new NetDyn;
    std::cout << "Running simulation " << i << "...\n";
    if (argc > 1)
    {
      configFile << argv[1];
      std::cout << "Loading Config File: " << configFile.str() << "\n";
      simulation->loadConfigFile(configFile.str(), i);
    }
    else
      simulation->loadConfigFile("config.cfg", i);
    i = simulation->simulationStart();
    seedsUsed << simulation->getOriginalRngSeed() << " ";
    // Bah. Need to avoid the memory hogs
    //        delete[] simulation;
  }
  std::cout << "List of simulation seeds:\n" << seedsUsed.str() << "\n";

  return 0;
}
