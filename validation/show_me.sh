for a in *.dat; do echo -e "\033[0;32m"$a"\033[0m";
	grep -A1 'result:' $a | grep 'result:' -v | cut -c 5-;
	echo -e "\033[0;31m"`grep -A1 'error:' $a | grep 'error:' -v | cut -c 5-`"\033[0m";
done;
