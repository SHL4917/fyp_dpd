filename=$(find . -mindepth 1 -maxdepth 1 -type d -printf '%f\n')

for item in $filename
do
	cd $item
	qsub combined_next.pbs
	cd ..
done
