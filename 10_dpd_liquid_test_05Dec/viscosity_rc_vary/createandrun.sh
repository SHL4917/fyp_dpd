for i in 0.6 0.75 0.8 0.9 1 1.05 1.1 1.2 1.3 2
do
	mkdir rc_$i
	cp *.pbs rc_$i
	cp in.* rc_$i	
	cd rc_$i
	rm *.txt
	rm *.pbs.*
	touch rc.txt
	echo $i >> rc.txt
		
	qsub *.pbs
	cd ..	
done
