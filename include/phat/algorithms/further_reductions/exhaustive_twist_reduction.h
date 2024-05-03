#pragma once

#include <phat/helpers/misc.h>
#include <phat/boundary_matrix.h>

namespace phat
{
    class exhaustive_twist_reduction
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

            for (index cur_dim = boundary_matrix.get_max_dim(); cur_dim >= 1; cur_dim--)
            {
                for (index cur_col = 0; cur_col < nr_columns; cur_col++)
                {
                    if (boundary_matrix.get_dim(cur_col) == cur_dim)
                    {
                        column final_col_value;
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
                            }
                            else
                            {
                                if (!is_pivot)
                                {
                                    lowest_one_lookup[lowest_one] = cur_col;
                                    boundary_matrix.clear(lowest_one);
                                    is_pivot = true;
                                }
                                final_col_value.push_back(lowest_one);
                                boundary_matrix.remove_max(cur_col);
                            }
                            lowest_one = boundary_matrix.get_max_index(cur_col);
                        }

                        column final_col_vector(final_col_value.rbegin(), final_col_value.rend());

                        boundary_matrix.set_col(cur_col, final_col_vector);
                        boundary_matrix.finalize(cur_col);
#if COUNT_OPS
                        count_final_ones += boundary_matrix.size(cur_col);
#endif
                    }
                }
            }
#if COUNT_OPS
	    phat::Reduction_statistics stats("Exhaustive-twist",count_ops,count_col_ops,count_final_ones,max_size);
	    stats.print();
#endif
        }
    };
}
