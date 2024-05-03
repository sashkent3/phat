#!/bin/sh

timeout=600

for combination in "--twist --bit_tree_pivot_column" "--swap_twist --full_pivot_column" "--exhaustive_compress --vector_vector" "--exhaustive_compress --bit_tree_pivot_column" "--retrospective --vector_vector"
do
echo INFO Running on ${combination} - primal
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii ${combination} --primal $1
echo INFO Running on ${combination} - dual
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii ${combination} --dual $1
done

