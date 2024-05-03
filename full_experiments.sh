#!/bin/sh

#remove timeout on bigger inputs to obtain the countings
timeout=600


#output is the fill-in, # of operations, and # of column operations of twist, swap, retro, exhaustive, and mix for the inputs
count_col_type=--vector_vector

echo "COUNTING OPERATIONS\n"
for input in ../inputs/alpha_*_010000_*.phat ../inputs/vr_random_50_16.phat ../inputs/vr_random_100_4.phat ../inputs/vr_senate.phat ../inputs/ls_nucleon_41x41x41_uint8.phat ../inputs/ls_fuel_64x64x64_uint8.phat ../inputs/ls_tooth_103x94x161_uint8.phat ../inputs/shuffled_50_*.phat ../inputs/shuffled_75_*.phat ../inputs/shuffled_100_*.phat
do
echo "INFO Create reference output"
./phat ${count_col_type} --ascii $input ref.out
echo "INFO Dualize input"
./convert --ascii --save-ascii --dualize $input $input.dual
echo "INFO Create reference output for dual"
./phat ${count_col_type} --ascii $input.dual ref_dual.out
for algorithm in --twist --swap_twist --retrospective --exhaustive_compress --mix_compress
do
rm -f bla
echo INFO Running on $algorithm
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout ${timeout} ./phat_with_count_ops --ascii $algorithm ${count_col_type} $input bla
diff -q ref.out bla
rm -f bla
echo INFO Running on $algorithm $column dualized
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout ${timeout} ./phat_with_count_ops --ascii $algorithm ${count_col_type} $input.dual bla
diff -q ref_dual.out bla
done
done

#output is the runtime over all data structures of twist, swap, retro, exhaustive, and mix for the inputs
echo "COMPARING DATASTRUCTURES\n"
for input in ../inputs/alpha_cube_040000_*.phat ../inputs/vr_celegans.phat ../inputs/ls_tooth_103x94x161_uint8.phat ../inputs/shuffled_75_*.phat
do
for algorithm in --twist --swap_twist --mix_compress --exhaustive_compress --retrospective
do
for column_type in --vector_list --vector_vector --vector_set --vector_heap --heap_pivot_column --sparse_pivot_column --full_pivot_column --bit_tree_pivot_column
do
echo INFO Running on $algorithm with ${column_type}
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii $algorithm ${column_type} --primal $input
echo INFO Running on $algorithm $column dualized with ${column_type}
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii $algorithm ${column_type} --dual $input
done
done


done

#output is the runtime over the best combination of data structures with twist, swap, retro, and exhaustive for the inputs
echo "COMPARING BEST COMBO\n"
for input in ../inputs/alpha_*_040000_*.phat ../inputs/alpha_*_080000_*.phat ../inputs/alpha_*_160000_*.phat ../inputs/vr_celegans.phat ../inputs/vr_house.phat ../inputs/vr_random_1000_8.phat ../inputs/ls_hydrogen_atom_128x128x128_uint8.phat ../inputs/ls_shockwave_64x64x512_uint8.phat ../inputs/ls_lobster_301x324x56_uint8.phat ../inputs/ls_mri_ventricles_256x256x124_uint8.phat ../inputs/ls_engine_256x256x128_uint8.phat ../inputs/ls_statue_leg_341x341x93_uint8.phat ../inputs/ls_bonsai_256x256x256_uint8.phat ../inputs/ls_skull_256x256x256_uint8.phat ../inputs/shuffled_1*_*.phat
do
for combination in "--twist --bit_tree_pivot_column" "--swap_twist --full_pivot_column" "--exhaustive_compress --vector_vector" "--exhaustive_compress --bit_tree_pivot_column" "--retrospective --vector_vector"
do
echo INFO Running on ${combination} - primal
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii ${combination} --primal $input
echo INFO Running on ${combination} - dual
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout $timeout ./benchmark --ascii ${combination} --dual $input
done
done





