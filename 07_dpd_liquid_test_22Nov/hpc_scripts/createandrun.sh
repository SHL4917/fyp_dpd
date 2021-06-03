for i in {1..7}
do
	mkdir $i
	cp diffusion.pbs $i
	cp in.diffusion_vacf $i	
	cd $i
	rm *.txt
	rm diffusion.pbs.*
		
	qsub diffusion.pbs
	cd ..	
done
