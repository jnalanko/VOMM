CXX = g++
STD = -std=c++11

.PHONY: bpr_to_dot score_string build_model build_model_optimized build_model_profile score_string_optimized tests maxreps_stats asd score_string_profile all profiling profiling tests

libraries= BD_BWT_index/lib/*.a sdsl-lite/build/lib/libsdsl.a sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a 
includes= -I BD_BWT_index/include -I sdsl-lite/include

all: tests score_string build_model
optimized: score_string_optimized build_model_optimized
profiling: score_string_profile build_model_profile

tests:
	$(CXX) $(STD) tests.cpp $(libraries) -o tests -Wall -Wno-sign-compare -Wextra $(includes) -g

bpr_to_dot:
	$(CXX) $(STD) bpr_to_dot.cpp -o bpr_to_dot -Wall -Wno-sign-compare -Wextra

build_model:
	$(CXX) $(STD) build_model.cpp $(libraries) -o build_model -Wall -Wno-sign-compare -Wextra $(includes) -g

build_model_optimized:
	$(CXX) $(STD) -O3 build_model.cpp $(libraries) -o build_model_optimized -Wall -Wno-sign-compare -Wextra $(includes) -g -march=native
	
build_model_profile:
	$(CXX) $(STD) build_model.cpp $(libraries) -o build_model_profile -Wall -Wno-sign-compare -Wextra $(includes) -O3 -g -pg
	
score_string:
	$(CXX) $(STD) score_string.cpp $(libraries) -o score_string -Wall -Wno-sign-compare -Wextra $(includes) -g

score_string_optimized:
	$(CXX) $(STD) -O3 score_string.cpp $(libraries) -o score_string_optimized -Wall -Wno-sign-compare -Wextra $(includes) -g -march=native
	
score_string_profile:
	$(CXX) $(STD) score_string.cpp $(libraries) -o score_string_profile -Wall -Wno-sign-compare -Wextra $(includes) -O3 -g -pg

maxreps_stats:
	$(CXX) $(STD) Maxreps_stats.cpp $(libraries) -O3 -o maxreps_stats -Wall -Wno-sign-compare -Wextra $(includes) -g
	
asd:
	$(CXX) $(STD) asd.cpp $(libraries) -o asd -Wall -Wno-sign-compare -Wextra $(includes) -g	
