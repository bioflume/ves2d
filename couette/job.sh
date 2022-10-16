#PBS -l nodes=1
#PBS -l walltime=168:00:00
#PBS -N p20_cntd
#PBS -m ae
#PBS -M rahimian@gatech.edu
#PBS -r y

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.02; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.02; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.03; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.03; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.04; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.04; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.05; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.05; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

#matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.06; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
#matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.06; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.07; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.07; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

#matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.08; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
#matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.08; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.09; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.09; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.1; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.1; reducedArea=.80; n=64; ts=1e-3; couette" > /dev/null 

matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.2; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.3; reducedArea=.65; n=64; ts=1e-3; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=1; volFrac=.4; reducedArea=.65; n=64; ts=5e-4; couette" > /dev/null &

# matlab -nodesktop -nodisplay -nosplash -r "viscCont=2; volFrac=.02; reducedArea=.65; n=128; ts=5e-5; RandSeed=30; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=2; volFrac=.02; reducedArea=.80; n=128; ts=5e-5; RandSeed=30; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=4; volFrac=.02; reducedArea=.65; n=128; ts=5e-5; RandSeed=30; couette" > /dev/null &
# matlab -nodesktop -nodisplay -nosplash -r "viscCont=4; volFrac=.02; reducedArea=.80; n=128; ts=5e-5; RandSeed=30; couette" > /dev/null &

wait 