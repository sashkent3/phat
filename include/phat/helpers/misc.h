/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

// STL includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <queue>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <iterator>

// VS2008 and below unfortunately do not support stdint.h
#if defined(_MSC_VER)&& _MSC_VER < 1600
    typedef __int8 int8_t;
    typedef unsigned __int8 uint8_t;
    typedef __int16 int16_t;
    typedef unsigned __int16 uint16_t;
    typedef __int32 int32_t;
    typedef unsigned __int32 uint32_t;
    typedef __int64 int64_t;
    typedef unsigned __int64 uint64_t;
#else
    #include <stdint.h>
#endif

// basic types. index can be changed to int32_t to save memory on small instances
namespace phat {
    typedef int64_t index;
    typedef int8_t dimension;
    typedef std::vector< index > column;
}

// OpenMP (proxy) functions
#if defined _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
	inline void omp_set_num_threads( int ) {};
    #include <time.h>
    #define omp_get_wtime() (float)clock() / (float)CLOCKS_PER_SEC
#endif

#include <phat/helpers/thread_local_storage.h>

#if COUNT_OPS

#include <locale>


// Number formatting taken from
// https://stackoverflow.com/questions/7276826/format-number-with-commas-in-c
class comma_numpunct : public std::numpunct<char>
{
  protected:
    virtual char do_thousands_sep() const
    {
        return ',';
    }

    virtual std::string do_grouping() const
    {
        return "\03";
    }
};

namespace phat {

  struct Reduction_statistics {

    std::string algname;
    long no_bitflips;
    long no_column_additions;
    long no_final_ones;
    long no_max_size;

    // Currently, only integer counts are supported
    std::vector<std::pair<std::string,long> > extra_fields; 

    void add_field(std::string key, long value) {
      extra_fields.push_back(std::make_pair(key,value));
    }

    Reduction_statistics(std::string _algname,
			 long _no_bitflips,
			 long _no_column_additions,
			 long _no_final_ones,
                         long _no_max_size) 
      : algname(_algname), no_bitflips(_no_bitflips),no_column_additions(_no_column_additions),no_final_ones(_no_final_ones),no_max_size(_no_max_size)
    {}

    void print() {
      std::locale comma_locale(std::locale(), new comma_numpunct());

      // tell cout to use our new locale.
      std::cout.imbue(comma_locale);
      
#if PRINT_LATEX
      
      std::cout << algname << " & " << no_bitflips << " & " << no_column_additions << " & " << no_final_ones << "\\\\" << std::endl;
#else
      std::cout << algname << std::endl;
      std::cout << std::left << std::setw(30) << "Nb bitflips" << ": \t" << std::right << std::setw(20) << no_bitflips << std::endl;
      std::cout << std::left << std::setw(30) << "Nb col ops" << ": \t" << std::right << std::setw(20) << no_column_additions << std::endl;
      std::cout << std::left << std::setw(30) << "Nb final 1s" << ": \t" << std::right << std::setw(20) << no_final_ones << std::endl;
      std::cout << std::left << std::setw(30) << "Max nb 1s" << ": \t" << std::right << std::setw(20) << no_max_size << std::endl;
      for(auto p : extra_fields) {
	std::cout << std::left << std::setw(30) << p.first << ": \t" << std::right << std::setw(20)<< p.second << std::endl;
      }
#endif
    }

  };

}
#endif

