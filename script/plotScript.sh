#!/bin/bash

#KOS 2016
#plot all data outputs from ksDft and generate laTeX document

#get all 2dMat files and plot them into eps - will be embedded in latex
MAINDIR=$PWD
echo "assemble gnuplot script!" | tee $MAINDIR/log.plotScript
echo "my current directory is" $PWD | tee -a tee $MAINDIR/log.plotScript

PLOTSCRIPT="plotscript.gnu"

MATRIXFILES=$(find -name "*.2dMat")

if [ -z "$MATRIXFILES" ]; then
  echo "no matrix files .2dMat found - exit!" | tee -a $MAINDIR/log.plotScript
  exit
else
 echo "found the following matrix files: " $MATRIXFILES  | tee -a $MAINDIR/log.plotScript
fi

POSTPROCDIR=$(find -name "POSTPROCESSING")
if [ -z "$POSTPROCDIR" ]; then
	echo "did not find postprocessing directory - create one!"  | tee -a $MAINDIR/log.plotScript
	mkdir POSTPROCESSING | tee -a $MAINDIR/log.plotScript
	POSTPROCDIR="POSTPROCESSING"
else
	echo "found the following postprocessing directory" $POSTPROCDIR  | tee -a $MAINDIR/log.plotScript
fi
DATADIR=$PWD"/POSTPROCESSING/DATA"
EPSDIR=$PWD"/POSTPROCESSING/EPS"
if [ ! -z "$DATADIR" ]; then
	  echo "did not find directory " $DATADIR " create one!"  | tee -a $MAINDIR/log.plotScript
	  mkdir ${DATADIR} | tee -a $MAINDIR/log.plotScript
else
	echo "found the following data directory" ${DATADIR}  | tee -a $MAINDIR/log.plotScript
fi
if [ ! -z "$EPSDIR" ]; then
	  echo "did not find directory " ${EPSDIR} " create one!"  | tee -a $MAINDIR/log.plotScript
	  mkdir ${EPSDIR} | tee -a $MAINDIR/log.plotScript
else
	echo "found the following eps directory" ${EPSDIR}  | tee -a $MAINDIR/log.plotScript
fi		

#now move MAT files into this folder 

for FILES in $MATRIXFILES 
do
  mv $FILES $POSTPROCDIR
done

#now change to postprocessing directory

cd $POSTPROCDIR

#now prepare the script

echo "we plot file",$PLOTFILENAME  | tee -a $MAINDIR/log.plotScript

echo "set xtic auto" | tee  $PLOTSCRIPT
echo "set ytic auto" | tee -a  $PLOTSCRIPT
echo "set title \"3d data\"" | tee -a  $PLOTSCRIPT
echo "set xlabel \"xlabel\"" | tee -a  $PLOTSCRIPT
echo "set ylabel \"ylabel\" " | tee -a  $PLOTSCRIPT
echo "set contour base" | tee -a $PLOTSCRIPT
echo "set cntrparam levels 7" | tee -a $PLOTSCRIPT
echo "unset surface" | tee -a $PLOTSCRIPT 
#echo "set palette rgbformula 33,13,10" | tee -a $PLOTSCRIPT
#echo "set size square" | tee -a $PLOTSCRIPT
echo "set view 0,0" | tee -a $PLOTSCRIPT
#echo "set palette maxcolors 12" | tee -a $PLOTSCRIPT
echo "splot inputFilename matrix with image" | tee -a $PLOTSCRIPT
echo "set terminal pdf " | tee -a $PLOTSCRIPT
echo "set output outputname " | tee -a  $PLOTSCRIPT
#echo "set output \"test.pdf\"" | tee -a $PLOTSCRIPT
echo "set contour base" | tee -a $PLOTSCRIPT
echo "set cntrparam levels 7" | tee -a $PLOTSCRIPT
echo "set contour base" | tee -a $PLOTSCRIPT
echo "set palette rgbformula 33,13,10" | tee -a $PLOTSCRIPT
#echo "set palette maxcolors 12" | tee -a $PLOTSCRIPT
echo "set view 0,0" | tee -a $PLOTSCRIPT
echo "splot inputFilename matrix with image" | tee -a $PLOTSCRIPT


#get all the files to plot and plot them all
for FILENAMEINLIST in $(find -name "*.2dMat")
do
  echo "plot file " $FILENAMEINLIST | tee -a $MAINDIR/log.plotScript
  FILENAME=$(basename $FILENAMEINLIST) 
  FILENAME_WO_ENDING=${FILENAME%%.*}
  OUTNAME=$FILENAME_WO_ENDING".pdf" 
  echo "The pdf outputfile is named " $OUTNAME | tee -a $MAINDIR/log.plotScript
  gnuplot -e "inputFilename='${FILENAMEINLIST}'; outputname='${OUTNAME}'" $PLOTSCRIPT  | tee -a $MAINDIR/log.plotScript
  echo "now move images to image directory"  | tee -a $MAINDIR/log.plotScript
  mv ${FILENAMEINLIST} $DATADIR
done

echo "move all pdf images to pdf directory"   | tee -a $MAINDIR/log.plotScript

for FILENAMEINLIST in $(find -name "*.pdf")
do
 mv $FILENAMEINLIST $EPSDIR
done

cd $MAINDIR

echo "now run PDFLATEX and recompile the latex document with the right figures!"   | tee -a $MAINDIR/log.plotScript

pdflatex latex.tex

echo "finished plotting! "   | tee -a $MAINDIR/log.plotScript
