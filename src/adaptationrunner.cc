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

#include "adaptationrunner.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>

adaptationRunner::adaptationRunner()
{
  initialize();
  baseMultiplier = 1.1;
  extrapolationMultiplier = 1;

  slopePositive = true;
  verbose = true;
}

adaptationRunner::adaptationRunner(double baseMult, double extrMult, bool slop, bool verb)
{
  initialize();
  baseMultiplier = baseMult;
  extrapolationMultiplier = extrMult;
  slopePositive = slop;
  verbose = verb;
}

void adaptationRunner::setLowerBound(double bound)
{
  yLowerBound = bound;
}

void adaptationRunner::setUpperBound(double bound)
{
  yUpperBound = bound;
}

void adaptationRunner::initialize()
{
  x1 = std::numeric_limits<double>::quiet_NaN();
  x2 = std::numeric_limits<double>::quiet_NaN();
  y1 = std::numeric_limits<double>::quiet_NaN();
  y2 = std::numeric_limits<double>::quiet_NaN();
  yLowerBound = -std::numeric_limits<float>::infinity();
  yUpperBound = std::numeric_limits<float>::infinity();
}

void adaptationRunner::addDataPair(double x, double y)
{
  // Add the data pair and shift the old values
  x1 = x2;
  y1 = y2;
  x2 = x;
  y2 = y;
  if (verbose)
    std::cout << "Adding adaptation point (" << x << "," << y << ").\n";
}

double adaptationRunner::getPrediction(double x)
{
  // First check
  if (std::isnan(x2) || std::isnan(y2))
  {
    std::cout << "No adaptation points available! Cannot predict anything.\n";
    exit(1);
  }
  // When only one point is present or there are infinities or the points are the same
  else if (std::isnan(x1) || std::isnan(y1) || std::isnan(x2) || std::isnan(y2) || std::isinf(x1) || std::isinf(x2) || x2 - x1 == 0.0 || y2 - y1 == 0.0)
  {
    if (verbose)
      std::cout << "Not enough datapoints available. Applying multiplier.\n";
    if (slopePositive)
    {
      if (x > x2)
        return y2 * baseMultiplier;
      else if (x < x2)
        return y2 / baseMultiplier;
      else
      {
        std::cout << "You are already at the target value!\n";
        return y2;
      }
    }
    else
    {
      if (x > x2)
        return y2 / baseMultiplier;
      else if (x < x2)
        return y2 * baseMultiplier;
      else
      {
        std::cout << "You are already at the target value!\n";
        return y1;
      }
    }
  }
  // When two points are present
  else
  {
    double interpolationRange = fabs(x2 - x1);
    double minX = std::min(x1, x2) - interpolationRange * extrapolationMultiplier;
    double maxX = std::max(x1, x2) + interpolationRange * extrapolationMultiplier;
    double newX = x;
    if (newX > maxX)
    {
      if (verbose)
        std::cout << "Desired value out of extrapolation bounds, going to the closest one.";
      newX = maxX;
    }
    if (newX < minX)
    {
      if (verbose)
        std::cout << "Desired value out of extrapolation bounds, going to the closest one.";
      newX = minX;
    }
    double newY = y1 + (y2 - y1) * (newX - x1) / (x2 - x1);
    if (newY > yUpperBound)
    {
      std::cout << "Predicted Y of " << newY << " above the upper bound. Something went wrong...";
      exit(1);
    }
    if (newY < yLowerBound)
    {
      std::cout << "Predicted Y of " << newY << " below the lower bound. Something went wrong...";
      exit(1);
    }
    if (verbose)
    {
      std::cout << "Predicted Y of " << newY << " for a target X of " << newX << "\n";
    }
    return newY;
  }
}
