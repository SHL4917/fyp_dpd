for i in {1..7}
do
	cd $i
	rm *.txt
	rm diffusion.pbs.*
	qsub diffusion.pbs
	cd ..	
done