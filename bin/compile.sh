cd ../build

cp ../src/Utils/kindSettings.f90 ./
cp ../src/Functions/funcs.f90 ./
cp ../src/Functions/derivatives.f90 ./
cp ../src/Functions/functionInterface.f90 ./
cp ../src/NumericalDerivatives/FiniteDifferences/centralDifferences.f90 ./
cp ../src/NumericalDerivatives/FiniteDifferences/forwardDifferences.f90 ./
cp ../src/NumericalDerivatives/FiniteDifferences/backwardDifferences.f90 ./
cp ../src/Test/test.f90 ./

gfortran -c kindSettings.f90
gfortran -c funcs.f90
gfortran -c derivatives.f90
gfortran -c functionInterface.f90
gfortran -c *Differences.f90

gfortran *.o test.f90
./a.out
rm a.out *.mod *.o
