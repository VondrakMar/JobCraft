#Thin uncomment if you want to take a filename as an argument


python ca

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 path/to/aimsRun.py /path/to/temp number_of_calcs_at_once number_of_lines_in_x_files"
    exit 1
fi
aimrRun=$1
temp=$2
at_one=$3
per_one_x=$4

for a in *xyz
do 
    mkdir ${a::-4}
    cp $a ${a::-4}
    cp $1 ${a::-4}
done

for a in */
do 
    echo "cd $a ; python aimsRun.py ${a::-1}.xyz" >> allCmds
done

split -$per_one_x allCmds

for a in x* 
do 
    cat $temp >> $a.sl ; 
    echo "parallel --delay 0.2 --joblog task.log --progress -j $at_one < $a" >> $a.sl 
done


#filename=$1
#if [ ! -f "$filename" ]; then
#    echo "File not found!"
#    exit 1
#fi
#
#

#
#for a in xa*; do cat temp >> $a.sl ; echo "parallel --delay 0.2 --joblog task.log --progress -j 9 < $a" >> $a.sl ; done
