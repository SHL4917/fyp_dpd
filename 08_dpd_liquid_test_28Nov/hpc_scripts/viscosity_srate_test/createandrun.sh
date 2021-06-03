for i in 0.05 0.10 0.15 0.18 0.25 0.3
do
	mkdir $i
	cp *.pbs $i
	cp in.* $i	
	cd $i
	rm *.txt
	rm *.pbs.*
	touch srate.txt
	echo $i >> srate.txt
		
	qsub *.pbs
	cd ..	
done
