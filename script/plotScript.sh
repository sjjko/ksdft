#!/bin/bash

#KOS 2016
MAINDIR=$PWD
echo "===========================================================================" | tee $MAINDIR/log.plotScript
echo "plotScript.sh: plot all data outputs from ksDft and generate laTeX document" | tee $MAINDIR/log.plotScript
echo "===========================================================================" | tee $MAINDIR/log.plotScript

if [ -z "$1" ]; then
  echo "Supply case name as first argument to this script - EXIT" | tee log.runsh
  echo "usage: bash ./plotScript.sh caseName" | tee log.runsh
  echo ""
  echo "the following case names are available:" | tee log.runsh
  CASELIST=$(echo case/*)$
  for CASE in $CASELIST
  do
    echo "		x " $(basename $CASE)
   done
  echo "===========================================================================" | tee $MAINDIR/log.plotScript
  exit
else
  CASENAME=$1
fi

#get all 2dMat files and plot them into eps - will be embedded in latex
echo "assemble gnuplot script!" | tee -a $MAINDIR/log.plotScript
echo "my current directory is" $PWD | tee -a $MAINDIR/log.plotScript


MATRIXFILES=$(find -maxdepth 1 -name "*.2dMat")

if [ -z "$MATRIXFILES" ]; then
  echo "no matrix files .2dMat found - exit!" | tee -a $MAINDIR/log.plotScript
  exit
else
 echo "found the following matrix files: " $MATRIXFILES  | tee -a $MAINDIR/log.plotScript
fi

POSTPROCDIR=$(find  -maxdepth 1 -name "POSTPROCESSING")
if [ ! -d "$POSTPROCDIR" ]; then
	echo "did not find postprocessing directory - create one!"  | tee -a $MAINDIR/log.plotScript
	mkdir POSTPROCESSING | tee -a $MAINDIR/log.plotScript
	POSTPROCDIR=$PWD"/POSTPROCESSING/"
else
	echo "found the following postprocessing directory" $POSTPROCDIR  | tee -a $MAINDIR/log.plotScript
fi

echo "now create a case directory for postprocessing data in " $CASEDIR | tee -a $MAINDIR/log.plotScript

POSTPROCDIR=$PWD"/POSTPROCESSING/"
CASEDIR=$POSTPROCDIR"/"$CASENAME"/"

if [ ! -d "$CASEDIR" ]; then
	echo "did not find case directory - create one!"  | tee -a $MAINDIR/log.plotScript
	mkdir $CASEDIR | tee -a $MAINDIR/log.plotScript
else
	echo "found a case directory" $CASEDIR  | tee -a $MAINDIR/log.plotScript
fi

DATADIR=$CASEDIR"/DATA/"
IMGDIR=$CASEDIR"/IMG/"

if [ ! -d "$DATADIR" ]; then
	  echo "did not find directory " $DATADIR " create one!"  | tee -a $MAINDIR/log.plotScript
	  mkdir ${DATADIR} | tee -a $MAINDIR/log.plotScript
else
	echo "found the following data directory" ${DATADIR}  | tee -a $MAINDIR/log.plotScript
fi
if [ ! -d "$IMGDIR" ]; then
	  echo "did not find directory " ${IMGDIR} " create one!"  | tee -a $MAINDIR/log.plotScript
	  mkdir ${IMGDIR} | tee -a $MAINDIR/log.plotScript
else
	echo "found the following eps directory" ${IMGDIR}  | tee -a $MAINDIR/log.plotScript
fi		

#now move MAT files into this folder 

for FILES in $MATRIXFILES #
do
  mv $FILES $CASEDIR
done

#now change to postprocessing directory

cd $CASEDIR

#now prepare the script

echo "we plot file",$PLOTFILENAME  | tee -a $MAINDIR/log.plotScript

PLOTSCRIPT="plotscript.gnu"

echo "set xtic auto" | tee  $PLOTSCRIPT
echo "set ytic auto" | tee -a  $PLOTSCRIPT
echo "set xlabel \"x [bohr]\"" | tee -a  $PLOTSCRIPT
echo "set ylabel \"y [bohr]\" " | tee -a  $PLOTSCRIPT
#echo "set contour base" | tee -a $PLOTSCRIPT
#echo "set cntrparam levels 7" | tee -a $PLOTSCRIPT
#echo "unset surface" | tee -a $PLOTSCRIPT 
#echo "set palette rgbformula 33,13,10" | tee -a $PLOTSCRIPT
#echo "set size square" | tee -a $PLOTSCRIPT
#echo "set view 0,0" | tee -a $PLOTSCRIPT
#echo "set palette maxcolors 12" | tee -a $PLOTSCRIPT
#echo "splot inputFilename matrix with image" | tee -a $PLOTSCRIPT
echo "set terminal pdf size 10 cm,10 cm" | tee -a $PLOTSCRIPT
echo "set output outputname " | tee -a  $PLOTSCRIPT
#echo "set output \"test.pdf\"" | tee -a $PLOTSCRIPT
echo "set contour base" | tee -a $PLOTSCRIPT
echo "set cntrparam levels 7" | tee -a $PLOTSCRIPT

echo "
set palette rgbformula 33,13,10
set cntrparam points 10
set cntrparam levels increment -6,-6,-24
set lmargin at screen 0.1
set rmargin at screen 0.85
set bmargin at screen 0.15
set tmargin at screen 0.85
unset border

" | tee -a $PLOTSCRIPT

echo "set contour base" | tee -a $PLOTSCRIPT
#echo "set palette rgbformula 33,13,10" | tee -a $PLOTSCRIPT
#  7,5,15   ... traditional pm3d (black-blue-red-yellow)
#  3,11,6   ... green-red-violet
#  23,28,3  ... ocean (green-blue-white); try also all other permutations
#  21,22,23 ... hot (black-red-yellow-white)
#  30,31,32 ... color printable on gray (black-blue-violet-yellow-white)
#  33,13,10 ... rainbow (blue-green-yellow-red)
#  34,35,36 ... AFM hot (black-red-yellow-white)
echo "set palette rgbformula 21,22,23" | tee -a $PLOTSCRIPT
#echo "set palette maxcolors 12" | tee -a $PLOTSCRIPT
#echo "set view 0,0" | tee -a $PLOTSCRIPT
#echo "set cntrparam bspline" | tee -a $PLOTSCRIPT
#echo "set table 'test.dat'
#splot inputFilename matrix
#unset table" | tee -a $PLOTSCRIPT
echo "unset key
set view map" | tee -a $PLOTSCRIPT

#echo "splot 'test.dat' with image" | tee -a $PLOTSCRIPT
echo "splot inputFilename matrix with image" | tee -a $PLOTSCRIPT


echo "set terminal pngcairo size 800,800 " | tee -a $PLOTSCRIPT
echo "set output outputnamepng " | tee -a $PLOTSCRIPT
echo "replot " | tee -a $PLOTSCRIPT

echo "set terminal pngcairo size 800,800 " | tee -a $PLOTSCRIPT
echo "set output outputnamepngWOAxis " | tee -a $PLOTSCRIPT
echo "unset ytics "| tee -a $PLOTSCRIPT
echo "unset xtics"| tee -a $PLOTSCRIPT
echo "unset border"| tee -a $PLOTSCRIPT
echo "unset colorbox"| tee -a $PLOTSCRIPT
echo "set xlabel \"\" " | tee -a  $PLOTSCRIPT
echo "set ylabel \"\" " | tee -a  $PLOTSCRIPT
echo "replot " | tee -a $PLOTSCRIPT


echo "we have the following files to apply gnuplot to: " $(find  -maxdepth 1 -name "*.2dMat")  | tee -a $MAINDIR/log.plotScript
for FILENAMEINLIST in $(find -maxdepth 1 -name "*.2dMat")
do
  echo "plot file " $FILENAMEINLIST | tee -a $MAINDIR/log.plotScript
  FILENAME=$(basename $FILENAMEINLIST) 
  FILENAME_WO_ENDING=${FILENAME%%.*}
  OUTNAME=$FILENAME_WO_ENDING".pdf" 
  OUTNAMEPNG=$FILENAME_WO_ENDING".png" 
  OUTNAMEPNGWOAXIS=$FILENAME_WO_ENDING"_woAxis.png" 
  echo "The pdf outputfile is named " $OUTNAME | tee -a $MAINDIR/log.plotScript
  gnuplot -e "inputFilename='${FILENAMEINLIST}'; outputname='${OUTNAME}'; outputnamepng='${OUTNAMEPNG}'; outputnamepngWOAxis='${OUTNAMEPNGWOAXIS}'" $PLOTSCRIPT  | tee -a $MAINDIR/log.plotScript
  echo "now move data to image directory"  | tee -a $MAINDIR/log.plotScript
  mv ${FILENAMEINLIST} $DATADIR
done

echo "move all pdf images to image directory " $IMGDIR  | tee -a $MAINDIR/log.plotScript

for FILENAMEINLIST in $(find  -maxdepth 1 -name "*.pdf")
do
 mv $FILENAMEINLIST $IMGDIR
done

for FILENAMEINLIST in $(find  -maxdepth 1 -name "*.png")
do
 mv $FILENAMEINLIST $IMGDIR
done

cd $MAINDIR

echo "now run PDFLATEX and recompile the latex document with the right figures!"   | tee -a $MAINDIR/log.plotScript

pdflatex latex.tex

echo "now cleanup and move latex files!"   | tee -a $MAINDIR/log.plotScript

latexFilename=$CASENAME"_doku.pdf"
echo "rename latex pdf to " $latexFilename  | tee -a $MAINDIR/log.plotScript

mv latex.pdf $latexFilename

LATEXDIR=$CASEDIR"doc"
echo "create new directory for documentation of the case in " $LATEXDIR | tee -a $MAINDIR/log.plotScript
mkdir $LATEXDIR

echo "move latex pdf to " $LATEXDIR | tee -a $MAINDIR/log.plotScript
mv $MAINDIR/$latexFilename $LATEXDIR

resultsFile="finalResult_"$CASENAME".dat"
echo "finally move final results file " $resultsFile " to results folder " $CASEDIR | tee -a $MAINDIR/log.plotScript
mv $MAINDIR/$resultsFile $CASEDIR

echo "now remove remnants of latex"  | tee -a $MAINDIR/log.plotScript
rm $MAINDIR/latex.*

echo "finished plotting! "   | tee -a $MAINDIR/log.plotScript
echo "===========================================================================" | tee -a $MAINDIR/log.plotScript
