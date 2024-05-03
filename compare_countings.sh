#!/bin/sh

#remove timeout on bigger inputs to obtain the countings
timeout=600

column_type=--vector_vector

echo "INFO Create reference output"
./phat ${column_type} --ascii $1 ref.out
echo "INFO Dualize input"
./convert --ascii --save-ascii --dualize $1 $1.dual
echo "INFO Create reference output for dual"
./phat ${column_type} --ascii $1.dual ref_dual.out
for algorithm in --twist --swap_twist --retrospective --exhaustive_compress --mix_compress
do
rm -f bla
echo INFO Running on $algorithm
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout ${timeout} ./phat_with_count_ops --ascii $algorithm ${column_type} $1 bla
diff -q ref.out bla
rm -f bla
echo INFO Running on $algorithm $column dualized
/usr/bin/time --format "INFO RESULTS %C\nINFO Time: %e\nINFO Memory: %M" timeout ${timeout} ./phat_with_count_ops --ascii $algorithm ${column_type} $1.dual bla
diff -q ref_dual.out bla
done

