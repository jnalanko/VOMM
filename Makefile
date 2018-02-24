.PHONY: bpr_to_dot score_string build_model build_model_optimized build_model_profile score_string_optimized tests maxreps_stats asd score_string_profile all profiling profiling

libraries= BD_BWT_index/lib/*.a sdsl-lite/build/lib/libsdsl.a sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a 
includes= -I BD_BWT_index/include -I sdsl-lite/include

all: tests score_string build_model
optimized: score_string_optimized build_model_optimized
profiling: score_string_profile build_model_profile

tests:
	g++ -std=c++11 tests.cpp $(libraries) -o tests -Wall -Wno-sign-compare -Wextra $(includes) -g

bpr_to_dot:
	g++ -std=c++11 bpr_to_dot.cpp -o bpr_to_dot -Wall -Wno-sign-compare -Wextra

build_model:
	g++ -std=c++11 build_model.cpp $(libraries) -o build_model -Wall -Wno-sign-compare -Wextra $(includes) -g

build_model_optimized:
	g++ -std=c++11 -O3 build_model.cpp $(libraries) -o build_model_optimized -Wall -Wno-sign-compare -Wextra $(includes) -g -march=native
	
build_model_profile:
	g++ -std=c++11 build_model.cpp $(libraries) -o build_model_profile -Wall -Wno-sign-compare -Wextra $(includes) -O3 -g -pg
	
score_string:
	g++ -std=c++11 score_string.cpp $(libraries) -o score_string -Wall -Wno-sign-compare -Wextra $(includes) -g

score_string_optimized:
	g++ -std=c++11 -O3 score_string.cpp $(libraries) -o score_string_optimized -Wall -Wno-sign-compare -Wextra $(includes) -g -march=native
	
score_string_profile:
	g++ -std=c++11 score_string.cpp $(libraries) -o score_string_profile -Wall -Wno-sign-compare -Wextra $(includes) -O3 -g -pg

maxreps_stats:
	g++ -std=c++11 Maxreps_stats.cpp $(libraries) -O3 -o maxreps_stats -Wall -Wno-sign-compare -Wextra $(includes) -g
	
asd:
	g++ -std=c++11 asd.cpp $(libraries) -o asd -Wall -Wno-sign-compare -Wextra $(includes) -g	
