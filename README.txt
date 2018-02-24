Compiling:

* Install sdsl-lite at the project root:

git clone https://github.com/simongog/sdsl-lite
cd sdsl-lite
sh install.sh
cd ..

* Install BD_BWT_index

cd BD_BWT_index
cmake -DCMAKE_BUILD_TYPE=Release .
make
cd ..

or

cd BD_BWT_index
cmake -DCMAKE_BUILD_TYPE=Debug . 
make
cd ..

* Compile rest with make

make tests


