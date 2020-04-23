
nearBragg :
	gcc -O -o bin/nearBragg nearBragg.c -lm -static
run :
	./bin/nearBragg -file dat/atoms.txt -lambda 1.0
	convert -depth 16 -type Grayscale -colorspace GRAY -endian LSB \
	-size 1024x1024+512 -negate -normalize \
	GRAY:intimage.img intimage.png
	mv floatimage.bin intimage.img intimage.png dat
