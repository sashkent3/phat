#!/bin/sh

timeout=600

for algorithm in --twist --swap_twist --mix_compress --exhaustive_compress --retrospective
do
for column_type in --vector_list --vector_vector --vector_set --vector_heap --heap_pivot_column --sparse_pivot_column --full_pivot_column --bit_tree_pivot_column 
do
echo INFO Running on $algorithm with ${column_type}
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii $algorithm ${column_type} --primal $1
echo INFO Running on $algorithm $column dualized with ${column_type}
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii $algorithm ${column_type} --dual $1
done
done

