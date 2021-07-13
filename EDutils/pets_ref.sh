file=$1
path=$(dirname $file)
tmp=tmp.txt
out=$path/out.txt
FoFc=$path/FoFc.txt

>$tmp
grep -A4 "^O[1-9].*ai" $file >> $tmp
grep -A4 "^C[1-9].*ai" $file >> $tmp
grep -A4 "^N[1-9].*ai" $file >> $tmp
grep -A4 "^H[1-9].*ai" $file >> $tmp


grep "^[OCNH]" $tmp | awk '{print $1}' >names.txt
grep "^ *3" $tmp | awk '{print $3, $4, $5}' > cols.txt

pr -mts' ' names.txt cols.txt > $out
rm names.txt cols.txt $tmp



echo 'h k l Fo Fc A B Fo-Fc sig(Fo) sqr(1/wt) sqr(wdFq) nref # SinThL ext sc' > $FoFc
# > $FoFc
sed -n '/Fc list after /,/Fc list of worst/p' $file \
| grep "^ *[0-9\-]" >> $FoFc
