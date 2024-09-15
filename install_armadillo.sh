# Shell scipt to automate build from source for armadillo library # 

PROJECT_ROOT=`pwd`

# Build OpenBLAS [Default Installation]
cd $PROJECT_ROOT/src/libraries/armadillo

echo "Building Armadillo"

cd src
sudo bash configure # Run Configuration File 
sudo cmake . # generate CMAKE Cache 
sudo cmake --build . # Build Armadillo from source

echo -e "{$GREEN} Armadillo build complete. Beginning Armadillo installation {$NC}"

sudo make install # Install armadillo 

GREEN='\033[0;32m'
NC='\033[0m' # No color

cd $PROJECT_ROOT # Return to Project Root 