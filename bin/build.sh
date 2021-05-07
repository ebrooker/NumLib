#!/bin/bash

NUMLIB_DIR=/data2/eb11d/misc/my_libs/github/NumLib/
cd $NUMLIB_DIR

if [ -d "$NUMLIB_DIR/build" ]
then
	echo "Build directory already exists, deleting it now..."
	rm -r build/
fi

echo "Creating build directory for CMake..."
mkdir build

echo "Entering the build directory..."
cd build

echo "Building NumLib project with CMake in BUILD MODE..."
cmake $NUMLIB_DIR/bin/

echo "Compiling NumLib project with make..."
make

if [ -f "$NUMLIB_DIR/build/numlib" ]
then
	echo "NumLib successfully built..."
	echo "Executing code..."
	./numlib
	echo ""
	echo "Have a nice day"
	echo ""
else
	echo "NumLib failed to build..."
	echo "What went wrong?"
	echo ""
fi