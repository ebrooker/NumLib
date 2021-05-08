
mkdir ../dev
cd ../dev

cp ../src/Utils/kindSettings.f90 ./
cp ../src/Functions/funcs.f90 ./
cp ../src/Functions/derivatives.f90 ./
cp ../src/Functions/functionInterface.f90 ./
cp ../src/Differentiation/RichardsonExtrapolation/richardsonExtrapolation.f90 ./
cp ../src/Differentiation/FiniteDifferences/finiteDifferences.f90 ./
cp ../src/Differentiation/FiniteDifferences/centralDifferences.f90 ./
cp ../src/Differentiation/FiniteDifferences/forwardDifferences.f90 ./
cp ../src/Differentiation/FiniteDifferences/backwardDifferences.f90 ./
cp ../src/Differentiation/differentiation.f90 ./
cp ../src/Test/proctorClass.f90 ./
cp ../src/Test/test.f90 ./

gfortran -c kindSettings.f90
gfortran -c funcs.f90
gfortran -c derivatives.f90
gfortran -c functionInterface.f90
gfortran -c richardsonExtrapolation.f90
gfortran -c forwardDifferences.f90
gfortran -c backwardDifferences.f90
gfortran -c centralDifferences.f90
gfortran -c finiteDifferences.f90
gfortran -c differentiation.f90
gfortran -c proctorClass.f90

gfortran *.o test.f90 -o numlib
./numlib
rm numlib *.mod *.o
