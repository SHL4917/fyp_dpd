for i in 1 1.001 1.01 1.015 1.02 1.03 1.05 1.1 1.15 1.2 1.3
do
	mkdir rc_$i
	cp *.pbs $i
	cp in.* $i	
	cd $i
	rm *.txt
	rm *.pbs.*
	touch rc.txt
	echo $i >> rc.txt
		
	qsub *.pbs
	cd ..	
done
