#!/bin/bash
MFILES="../../src/*.m"

#clear old mfiles.tex if it exists
if [ -f mfiles.tex ]
then
	rm mfiles.tex
fi

#Add code listing commands
for f in $MFILES
do
	mname=`basename $f`
	mname=`echo $mname | sed 's/_/\\\_/g'`
	echo "\subsection{$mname}" >> mfiles.tex
	echo "\lstinputlisting{$f}" >> mfiles.tex
done
