#!/bin/bash

#####################################
#           VARIABLES               #
#####################################

#Function to be integrated
MY_CFUN=$1
#Cuba c++ code to integrate the function
CUBA_CFUN=cuba.c
CUBA_METHOD=1 #1:Vegas 2:Suave 3:Divone 4:Cuhre

#Cuba Libraries
case "$OSTYPE" in
    darwin*)  echo "OSX" 
	      CUBA_FOLDER=/usr/local/Cellar/cuba/4.2
	      MY_LIB=$CUBA_FOLDER/lib/
	      MY_INC=$CUBA_FOLDER/include/	
	      ;;
    linux*)   echo "LINUX" 
	      #CUBA_FOLDER=:~/Programs/Cuba-4.2:~/Cuba-4.2
	      CUBA_FOLDER=~/Programs/Cuba-4.2
	      MY_LIB=$CUBA_FOLDER/
	      MY_INC=$CUBA_FOLDER/
	      ;;
    
    *)        echo "unknown: $OSTYPE" ; exit 2;;
esac




#####################################
#             COMPILE               #
#####################################
echo "g++ -std=c++11 -O3 $MY_CFUN.cpp -c -I$MY_INC"
# DEBUG
#export CUBACORES=1
if [ $# -eq 2 ]; then 
    g++ -std=c++11 -O3 -D ROTATION_ANGLE=$2 $MY_CFUN.cpp -c -I$MY_INC
else
    g++ -std=c++11 -O3 $MY_CFUN.cpp -c -I$MY_INC
fi
g++ -std=c++11 -O3 -fPIC DCD.cpp -c -I$MY_INC 

#Compile Cuba
g++ -std=c++11 -O4 $CUBA_CFUN -D CUBA_WHAT=$CUBA_METHOD -I$MY_INC -L$MY_LIB -lcuba $MY_CFUN.o DCD.o
