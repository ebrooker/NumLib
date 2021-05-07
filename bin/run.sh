
NUMLIB_DIR=/data2/eb11d/misc/my_libs/github/NumLib/

cd $NUMLIB_DIR/build

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