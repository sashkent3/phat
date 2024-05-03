#pragma once

#include <phat/helpers/misc.h>
#include <phat/boundary_matrix.h>
#include <stack>

namespace phat {
    class lazy_exhaustive_reduction {
    private:
      

      std::vector< index > positivePair;
      std::vector< index > negativePair;
      

#if COUNT_OPS
      long count_ops;
      long count_col_ops;
      long count_final_ones;
      long count_forward_col_ops;
      long count_forward_bitflips;
      long count_backward_col_ops;
      long count_backward_bitflips;

      std::vector<bool> used_for_addition;
#endif

    public:

      template< typename Representation >
	void reduce_exhaustively(boundary_matrix< Representation >& boundary_matrix,
				 index i) {
	
	column result;
        while(!boundary_matrix.is_empty(i)) {
	  index next = boundary_matrix.get_max_index(i);
	  index pair = positivePair[next];
	  if(pair==-1 || pair==i) {
	    result.push_back(next);
	    boundary_matrix.remove_max(i);
	  } else {
	    reduce_exhaustively(boundary_matrix,pair);
	    boundary_matrix.add_to(pair,i);
#if COUNT_OPS
	    count_ops += boundary_matrix.size(pair);
	    count_col_ops += 1;
	    if(pair<i) {
	      count_forward_col_ops++;
	      count_forward_bitflips+=boundary_matrix.size(pair);
	    } else {
	      count_backward_col_ops++;
	      count_backward_bitflips+=boundary_matrix.size(pair);
	    } 
	    used_for_addition[pair]=true;
#endif
	  }
	}
	std::reverse(result.begin(),result.end());
	boundary_matrix.set_col(i,result);
      }
    
      

      template< typename Representation >
        void operator () ( boundary_matrix< Representation >& boundary_matrix ) {
        

            const index nr_columns = boundary_matrix.get_num_cols();

#if COUNT_OPS
            count_ops = 0;
            count_col_ops = 0;
            count_final_ones = 0;
	    count_forward_col_ops=0;
	    count_forward_bitflips=0;
	    count_backward_col_ops=0;
	    count_backward_bitflips=0;

	    used_for_addition.resize(nr_columns,false);
#endif

	    

            positivePair.resize( nr_columns, -1 );
            negativePair.resize( nr_columns, -1 );
            
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


                index lowest_one = boundary_matrix.get_max_index(j);
                while (lowest_one != -1 && positivePair[lowest_one] != -1)
                {
		  index col_to_add=positivePair[lowest_one];
		  reduce_exhaustively(boundary_matrix,col_to_add);
#if COUNT_OPS
                    count_ops += boundary_matrix.size(col_to_add);
                    count_col_ops += 1;
		    count_forward_col_ops++;
		    count_forward_bitflips+=boundary_matrix.size(col_to_add);
		    used_for_addition[col_to_add]=true;
#endif
                    boundary_matrix.add_to(col_to_add, j);
                    lowest_one = boundary_matrix.get_max_index(j);
                }
                if (lowest_one != -1)
                {
		  index pair = boundary_matrix.get_max_index(j);
		  positivePair[pair] = j;
		  negativePair[j] = pair;
		}
            }

#if COUNT_OPS
	    long nr_used_for_addition=0;
            for( index j = 0; j < nr_columns; j++ ) {
                count_final_ones += boundary_matrix.size(j);
		if(used_for_addition[j]) {
		  nr_used_for_addition++;
		}
            }
#endif
            
#if COUNT_OPS
            long max_size = 0;
	    phat::Reduction_statistics stats("Lazy-retrospective",count_ops,count_col_ops,count_final_ones,max_size);
	    stats.add_field("Forward col ops",count_forward_col_ops);
	    stats.add_field("Forward bitflips",count_forward_bitflips);
	    stats.add_field("Backward col ops",count_backward_col_ops);
	    stats.add_field("Backward bitflips",count_backward_bitflips);
	    stats.add_field("Number of columns",nr_columns);
	    stats.add_field("Used for addition",nr_used_for_addition);
	    stats.print();

#endif

	}	    

    };
}
