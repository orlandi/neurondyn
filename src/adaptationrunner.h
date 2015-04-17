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

class adaptationRunner
{
  public:
  adaptationRunner();
  adaptationRunner(double baseMult, double extrMult, bool slop, bool verb = true);
  void initialize();
  void setLowerBound(double bound);
  void setUpperBound(double bound);
  void addDataPair(double x, double y);
  double getPrediction(double x);

  private:
  double x1, x2, y1, y2, yLowerBound, yUpperBound;
  double baseMultiplier, extrapolationMultiplier;
  bool slopePositive, verbose;
};