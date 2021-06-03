for i in 10 15 20 25 30 35
do
	cd e$i		
	qsub polymer.pbs

	cd ..	
done
