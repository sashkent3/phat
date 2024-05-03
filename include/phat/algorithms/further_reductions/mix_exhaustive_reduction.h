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
    class mix_exhaustive_reduction
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

#if COUNT_OPS
            long max_size = 0;
            std::vector<index> sizes(nr_columns, 0);
            for (index cur_col = 0; cur_col < nr_columns; cur_col++)
            {
                sizes[cur_col] = boundary_matrix.size(cur_col);
                max_size += sizes[cur_col];
            }

#endif

            for (index cur_col = 0; cur_col < nr_columns; cur_col++)
            {
                int final_col_value_nb = 0; // Current nb of 1s before the "lowest_one" index
                index sparsiest_index;      // Index of the sparsiest column observed so far
                index sparsiest_num;        // Nb of 1s in the sparsiest column
                column full_pivot_column;   // Column value that the std reduction would have returned (used as checkpoint)
                std::list<index> add_ops;   // List of add operations performed during exhaustive reduction

                bool is_pivot = false;
                index lowest_one = boundary_matrix.get_max_index(cur_col);
                while (lowest_one != -1)
                {
                    if (lowest_one_lookup[lowest_one] != -1)
                    {
                        #if COUNT_OPS
                        count_ops += boundary_matrix.size(lowest_one_lookup[lowest_one]);
                        count_col_ops += 1;
                        #endif

                        boundary_matrix.add_to(lowest_one_lookup[lowest_one], cur_col);

                        if (is_pivot)
                        { // Check if current column is sparser than pervious sparsiest one
                            add_ops.insert(add_ops.end(), lowest_one);
                            int cur_nb_ones = final_col_value_nb + boundary_matrix.size(cur_col);
                            if (cur_nb_ones < sparsiest_num)
                            {
                                sparsiest_index = lowest_one;
                                sparsiest_num = cur_nb_ones;
                            }
                        }
                    }
                    else
                    {
                        if (!is_pivot)
                        {
                            lowest_one_lookup[lowest_one] = cur_col;
                            is_pivot = true;

                            sparsiest_index = lowest_one;
                            sparsiest_num = boundary_matrix.size(cur_col);
                            boundary_matrix.get_col(cur_col, full_pivot_column);
                        }
                        final_col_value_nb += 1;
                        boundary_matrix.remove_max(cur_col);
                    }
                    lowest_one = boundary_matrix.get_max_index(cur_col);
                }
                // Reconstruct the sparsiest column
                if (is_pivot)
                {
                    boundary_matrix.set_col(cur_col, full_pivot_column);
                    for (index op_index : add_ops)
                    {
                        if (op_index < sparsiest_index)
                        {
                            break;
                        }
                        boundary_matrix.add_to(lowest_one_lookup[op_index], cur_col);
                    }
                }
                boundary_matrix.finalize(cur_col);
                #if COUNT_OPS
                count_final_ones += boundary_matrix.size(cur_col);
                #endif
            }
            #if COUNT_OPS
	    phat::Reduction_statistics stats("Mix-exhaustive",count_ops,count_col_ops,count_final_ones,max_size);
	    stats.print();
            #endif
        }
    };
}
