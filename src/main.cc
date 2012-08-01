/*
 * Copyright (c) 2009-2010 Javier G. Orlandi <orlandi@dherkova.com>
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

#include <QApplication>
#include <iostream>
#include <sstream>
#include "main.h"
#include "netdyn.h"

int main(int argc, char *argv[])
{
    int i = 1;
    std::stringstream seedsUsed;
    while(i > 0)
    {
        simulation = new NetDyn;
        std::cout << "Running simulation " << i << "...\n";
        simulation->loadConfigFile("config.cfg", i);
        i = simulation->simulationStart();
        seedsUsed << simulation->getOriginalRngSeed() << " ";
        // Bah. Need to avoid the memory hogs
//        delete[] simulation;
    }
    std::cout << "List of simulation seeds:\n" << seedsUsed.str() << "\n";

    return 0;
}

