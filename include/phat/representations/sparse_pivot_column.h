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

#include <phat/helpers/misc.h>

namespace phat {
    class sparse_column {

    protected:
        std::set< index > data;

        void add_index( const index idx ) {
            std::pair< std::set< index >::iterator, bool > result = data.insert( idx );
            if( result.second == false )
                data.erase( result.first );
        }

    public:
        void init( const index total_size ) {
            data.clear(); 
        }
	
	
	template<typename InputIterator>
	  void add_col(InputIterator begin, InputIterator end) {
	  for(InputIterator it=begin;it!=end;it++) {
	      add_index(*it);
	  }
	}
	

        void add_col( const column& col ) {
            for( index idx = 0; idx < (index) col.size(); idx++ )
                add_index( col[ idx ] );
        }

        index get_max_index() {
            return data.empty() ? -1 : *data.rbegin();
        }

        void get_col_and_clear( column& col ) {
            col.assign( data.begin(), data.end() );
            data.clear();
        }

        bool is_empty() {
            return data.empty();
        }

	void clear() {
	    data.clear();
	}

	index size() {
   	    return static_cast<index>(data.size());
	}
		
	void remove_max() {
            add_index( get_max_index() );
        }

        void set_col( const column& col  ) {
            clear();
            add_col( col );
        }

        void get_col( column& col  ) {
            get_col_and_clear( col );
            add_col( col );
        }
    };

}
