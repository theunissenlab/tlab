sex=$1
bird=$2
guess=$3
name=$sex'bird'$2
for birdwav in ../data/$name'call'*.wav;
do
	python fitCallFull.py $birdwav $guess > ../python_output/$birdwav'.fit'
done
