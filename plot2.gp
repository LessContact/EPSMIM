nx=801
ny=801
filename="cmake-build-release/double02000.dat"
set terminal png size 700,600
set output filename.".png"
set xrange[-1:nx]
set yrange[-1:ny]
set zrange[-4e-5:4e-5]
set cbrange[-4e-5:4e-5]
set palette gray
set title filename
plot filename binary array=(ny,nx) format="%lf" with image