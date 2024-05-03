/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Talha Bin Masood, Barbara Giunti, Guillaume Houry, Michael Kerber, Abhishek Rathod

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
#include <phat/boundary_matrix.h>

namespace phat
{
    class compress_reduction
    {
    public:
        template <typename Representation>
        void operator()(boundary_matrix<Representation> &boundary_matrix)
        {
#if COUNT_OPS
            long count_ops = 0;
            long count_col_ops = 0;
            long count_final_ones = 0;
#endif



            const index nr_columns = boundary_matrix.get_num_cols();
            std::vector<index> lowest_one_lookup(nr_columns, -1);

	    std::vector<bool> negative(nr_columns,false);

            for (index cur_col = 0; cur_col < nr_columns; cur_col++)
            {

	        // Compress
	        column faces;
                boundary_matrix.get_col(cur_col, faces);
                column boundary;
                for( index face : faces ) {
                    if(!negative[face]) {
                        boundary.push_back(face);
                    }
                }
                boundary_matrix.set_col(cur_col, boundary);


                index lowest_one = boundary_matrix.get_max_index(cur_col);
                while (lowest_one != -1 && lowest_one_lookup[lowest_one] != -1)
                {
#if COUNT_OPS
                    count_ops += boundary_matrix.size(lowest_one_lookup[lowest_one]);
                    count_col_ops += 1;
#endif
                    boundary_matrix.add_to(lowest_one_lookup[lowest_one], cur_col);
                    lowest_one = boundary_matrix.get_max_index(cur_col);
                }
                if (lowest_one != -1)
                {
                    lowest_one_lookup[lowest_one] = cur_col;
		    negative[cur_col]=true;
                }
                boundary_matrix.finalize(cur_col);
#if COUNT_OPS
                count_final_ones += boundary_matrix.size(cur_col);
#endif
            }
#if COUNT_OP
	    phat::Reduction_statistics stats("Compress",count_ops,count_col_ops,count_final_ones);
	    stats.print();
#endif
        }
    };
}
