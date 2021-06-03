for i in {4..14}
do
	cd e$i		
	qsub polymer.pbs

	cd ..	
done
