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
#include <stack>

namespace phat {
    class greedy_retrospective_reduction {
    public:
        template< typename Representation >
        void operator () ( boundary_matrix< Representation >& boundary_matrix ) {

#if COUNT_OPS
            long count_ops = 0;
            long count_col_ops = 0;
            long count_final_ones = 0;
#endif

            const index nr_columns = boundary_matrix.get_num_cols();
            std::vector< index > positivePair( nr_columns, -1 );
            std::vector< index > negativePair( nr_columns, -1 );

#if COUNT_OPS
            long max_size = 0;
            std::vector<index> sizes(nr_columns, 0);
            for (index cur_col = 0; cur_col < nr_columns; cur_col++)
            {
                sizes[cur_col] = boundary_matrix.size(cur_col);
                max_size += sizes[cur_col];
            }

#endif

            std::vector< std::vector< index > > boundary_transpose(nr_columns);
            for( index j = 0; j < nr_columns; j++ ) {
                if (boundary_matrix.is_empty(j)) {
                    continue;
                }
                column faces;
                boundary_matrix.get_col(j, faces);
                for( index face : faces ) {
                    boundary_transpose[face].push_back(j);
                }
            }

            for( index j = 0; j < nr_columns; j++ ) {
                if (boundary_matrix.is_empty(j)) {
                    continue;
                }

                // Compress
                column faces;
                boundary_matrix.get_col(j, faces);
                column boundary;
                for( index face : faces ) {
                    if(negativePair[face] == -1) {
                        boundary.push_back(face);
                    }
                }
                boundary_matrix.set_col(j, boundary);

                for( index face : boundary) {
                    index pair = positivePair[face];
                    if( pair != -1 && pair != j) {
                        boundary_matrix.add_to(pair, j);
#if COUNT_OPS
                        count_ops += boundary_matrix.size(pair);
                        count_col_ops += 1;
#endif
                    }
                }
                column new_col;
                boundary_matrix.get_col(j, new_col);
                update_transpose(faces, new_col, j, boundary_transpose);

                if (!boundary_matrix.is_empty(j)) {
                    index pair = boundary_matrix.get_max_index(j);
                    positivePair[pair] = j;
                    negativePair[j] = pair;

                    std::vector< index > row_copy(boundary_transpose[pair]);

                    for (index col: row_copy) {
                        if (col != j){
                            column prev_col;
                            boundary_matrix.get_col(col, prev_col);
                            boundary_matrix.add_to(j, col);
#if COUNT_OPS
                            count_ops += boundary_matrix.size(j);
                            count_col_ops += 1;
#endif
                            column new_col;
                            boundary_matrix.get_col(col, new_col);
                            update_transpose(prev_col, new_col, col, boundary_transpose);
                        }
                    }
                }
            }

#if COUNT_OPS
            for( index j = 0; j < nr_columns; j++ ) {
                count_final_ones += boundary_matrix.size(j);
            }

	    phat::Reduction_statistics stats("Greedy-retrospective",count_ops,count_col_ops,count_final_ones,max_size);
	    stats.print();
#endif
        }

    private:
        void update_transpose(column &prev, column &curr, index col_idx, std::vector< std::vector< index > > &boundary_transpose) {
            std::sort(prev.begin(), prev.end());
            std::sort(curr.begin(), curr.end());
            auto prev_col_it = prev.begin();
            auto curr_col_it = curr.begin();
            while (prev_col_it != prev.end() && curr_col_it != curr.end()){
                if (*prev_col_it == *curr_col_it){
                    prev_col_it++;
                    curr_col_it++;
                } else if (*prev_col_it < *curr_col_it) {
                    index face = *prev_col_it;
                    auto pr = std::equal_range(std::begin(boundary_transpose[face]), std::end(boundary_transpose[face]), col_idx);
                    boundary_transpose[face].erase(pr.first, pr.second);
                    prev_col_it++;
                } else {
                    index face = *curr_col_it;
                    boundary_transpose[face].insert( std::upper_bound( boundary_transpose[face].begin(), boundary_transpose[face].end(), col_idx ), col_idx );
                    curr_col_it++;
                }
            }
            while (prev_col_it != prev.end()) {
                index face = *prev_col_it;
                auto pr = std::equal_range(std::begin(boundary_transpose[face]), std::end(boundary_transpose[face]), col_idx);
                boundary_transpose[face].erase(pr.first, pr.second);
                prev_col_it++;
            }
            while (curr_col_it != curr.end()) {
                index face = *curr_col_it;
                boundary_transpose[face].insert( std::upper_bound( boundary_transpose[face].begin(), boundary_transpose[face].end(), col_idx ), col_idx );
                curr_col_it++;
            }
        }
    };
}
