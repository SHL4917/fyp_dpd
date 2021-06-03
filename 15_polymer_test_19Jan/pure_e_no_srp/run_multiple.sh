for i in {5..17}
do
	cd e$i		
	qsub polymer.pbs

	cd ..	
done
