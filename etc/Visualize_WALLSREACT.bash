#!/bin/bash
#author: Ron Ballouz
#Date: Winter/Spring 2015-2016
#Description:
#Script to Generate movie for WALLS_REACT simulation
#note: if this script crashes or exits before completing make sure to execute the following two commands before restarting:
#     cp ssdraw_def.par ssdraw.par
#     cp walls_def.dat walls.dat

frame_rate=30

if [ ! -f walls.dat ]; then
  echo "walls.dat file not found!"
  exit 1
fi

if [ ! -f ssdraw.par ]; then
  echo "ssdraw.par file not found!"
  exit 1
fi

if [ ! -f povray.inc ]; then
  echo "povray.inc file not found!"
  exit 1
fi

if ls ss.0*r 1> /dev/null 2>&1; then
    touch /dev/null
else
  echo "reduced output files not found!"
  exit 1
fi

cp walls.dat walls_def.dat
cp ssdraw.par ssdraw_def.par

################Python Step
echo "Python Step - Creating wall and ssdraw files"
python movewallsFileCreator.py
wait

##############Drawing Individual Frames
echo "Making .PNG's"
red_files=`ls ss.*.r`

for i in $red_files
do
  frame=`echo $i | cut -d. -f2`
  cp ssdraw$frame.par ssdraw.par
  ssdraw $i&
  wait
  pov ${i%*.ss}.pov.gz >& pov_temp& #replace .ss with .pov in string i
  wait
  rm pov_temp
done

#Cleanup
rm *.pov.gz walls0*.dat ssdraw0*.par

wait

##############Making the Movie
echo "Making Movie"

pngs=`\ls -1d ss.*.r.png`
i=1
for png in $pngs
do
  ln -s $png ffmpeg`printf %012d $i`.png
  i=$((i+1))
done
rm movie.mp4
ffmpeg -i ffmpeg%012d.png -r $frame_rate -pix_fmt yuv420p movie.mp4 << EOF

EOF

####Cleaning up Again
rm -f ffmpeg*.png >& /dev/null

cp ssdraw_def.par ssdraw.par
cp walls_def.dat walls.dat

rm ssdraw_def.par walls_def.dat ss.0*png
