path=$1
name=$2
pathname=$path/$name
out=$path"/dat"

if [ ! -d $out ];then echo "creating $out"; mkdir $out;fi

echo "...pts..."
sed -n '/lambda/,/List/p' $pathname.pts \
| grep "lambda\|omega\|aperpixel" > $out/pts.txt #\
# | awk '{print $2;}'


sed -n '/imagelist/,/endimagelist/p' $pathname.pts \
| sed '1d; $d' |  sed 's/\( *\) /,/g' \
| sed 's/^,//'> $out/iml.txt



echo "...rpl xyz cor hkl cenloc..."
types="rpl xyz cor hkl cenloc"
for type in $types; do
    sed 's/\( *\) /,/g' $pathname.$type | sed 's/^,//'  > $out/$type.txt
done


echo "...cif..."
sed -n '/_diffrn_zone_axis_scale/,/loop_/p'  $pathname.cif_pets \
| sed '1d; $d' |    sed 's/\( *\) /,/g' \
| sed 's/^,//'  > $out/cif.txt

sed -n '/_refln_zone_axis_id/,//p'  $pathname.cif_pets \
| sed '1d' |    sed 's/\( *\) /,/g' \
| sed 's/^,//'  > $out/HKL.txt



echo "...orientation matrix..."
grep _diffrn_orient_matrix_UB $pathname.cif_pets | awk '{print $2}' > $out/UB.txt
echo "...lattice_parameters..."
grep _cell_ $pathname.cif_pets | awk '{print $2}' > $out/cell.txt


dyn=$pathname"_dyn.cif_pets"
echo "...dyn..."
sed -n '/_diffrn_zone_axis_scale/,/loop_/p'  $dyn \
| sed '1d; $d' |    sed 's/\( *\) /,/g' \
| sed 's/^,//'  > $out/dyn.txt

sed -n '/_refln_zone_axis_id/,//p'  $dyn \
| sed '1d' |    sed 's/\( *\) /,/g' \
| sed 's/^,//'  > $out/HKL_dyn.txt


echo "...orientation matrix..."
grep _diffrn_orient_matrix_UB $dyn | awk '{print $2}' > $out/UB_dyn.txt
echo "...lattice_parameters..."
grep _cell_ $dyn | awk '{print $2}' > $out/cell_dyn.txt


echo "...fix hkl S line..."
sed -i 's/,S//'  $out/hkl.txt
