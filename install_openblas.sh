# Shell scipt to automate build from source for openblas library # 

PROJECT_ROOT=`pwd`

# Build OpenBLAS [Default Installation]
cd $PROJECT_ROOT/src/libraries/openblas

if [! -d "build" ]; then 
    echo "Building OpenBLAS"
fi

if [ -d "build" ]; then 
    rm -rf build # Clear Previous Installation
    echo "Rebuilding OpenBLAS"
fi 

mkdir build 
cd build 
sudo cmake ../src # Generate CMAKE Cache 
sudo cmake --build . # Build openblas

GREEN='\033[0;32m'
NC='\033[0m' # No color

echo -e "{$GREEN} OpenBLAS build complete. Beginning OpenBLAS installation {$NC}"

sudo make PREFIX=. install # Install openblas

cd $PROJECT_ROOT # Return to project root 