for ((bird=1;bird<11;bird++));
do
	for((call=1;call<16;call++));
	do
		python fit_bga.py m $bird $call > output-m-$bird-$call.txt
		python fit_bga.py m $bird $call
	done
done

