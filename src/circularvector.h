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

#pragma once
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <vector>
#include <cmath>

// I cannot believe this program didn't use the stl containers...
class circularVector
{
  public:
  circularVector(int s = 1000)
  {
    size = s;
    idx = 0;
    vec = std::vector<double>(size);
    sumX = 0.;
    sumXsquared = 0.;
    oldValue = 0.;
    vectorFull = false;
  }
  void increaseCurrentBin()
  {
    vec[idx]++;
  }
  void addToCurrentBin(double value)
  {
    vec[idx] += value;
  }
  double getCurrentValue()
  {
    return vec[idx];
  }
  void moveToNextBin()
  {
    //idx++;
    // Updated the sums (remove the old value and add the new)
    sumX += vec[idx] - oldValue;
    sumXsquared += vec[idx] * vec[idx] - oldValue * oldValue;
    idx++;
    // No need to keep track of the global index
    if (idx == size)
    {
      vectorFull = true;
      idx = 0;
    }
    oldValue = vec[idx];
    vec[idx] = 0;
    //idx % size
  }
  bool isVectorFull()
  {
    return vectorFull;
  }
  double getMean()
  {
    return sumX / double(size);
  }
  double getStdDeviation()
  {
    return sqrt(double(size) * sumXsquared - sumX * sumX) / double(size);
  }
  void clear()
  {
    std::fill(vec.begin(), vec.end(), 0);
    idx = 0;
    sumX = 0.;
    sumXsquared = 0.;
    oldValue = 0.;
    vectorFull = false;
  }

  private:
  bool vectorFull;
  int idx, size;
  std::vector<double> vec;
  double oldValue;
  double sumX, sumXsquared;
};