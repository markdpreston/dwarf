CXXFLAGS  = -std=c++0x -pedantic -ansi -O3 -Wall -I. -I/usr/include/gsl $(shell /usr/bin/R CMD config --cppflags) $(shell echo 'Rcpp:::CxxFlags()' | /usr/bin/R --vanilla --slave) $(shell echo 'RInside:::CxxFlags()' | /usr/bin/R --vanilla --slave)
VPATH     = .
SRC       = $(wildcard *.cpp)
OBJ       = $(patsubst %.cpp,.obj/%.o,$(SRC))
LIB       = -lm -l:/usr/lib/libgsl.a -l:/usr/lib/libgslcblas.a -l:/usr/lib/libf2c.a -l:/usr/lib/libboost_system.a -L/usr/lib/R/lib -lR -l:/usr/local/lib/R/site-library/RInside/lib/libRInside.a -l:/usr/local/lib/R/site-library/Rcpp/lib/libRcpp.a

all : $(OBJ) .obj/sadmvnt.o
	@echo "Linking..."
	@$(CXX) $(CXXFLAGS) -o ../dwarf $(OBJ) .obj/sadmvnt.o $(LIB)

.obj/sadmvnt.o : sadmvnt.f
	@echo $<
	@f77 -c sadmvnt.f -o .obj/sadmvnt.o 2> /dev/null

.obj/%.o : %.cpp
	@echo $<
	@mkdir -p .obj
	$(CXX) $(CXXFLAGS) -c $< -o $(patsubst %.cpp,.obj/%.o,$<)

clean :
	@rm -f ../dwarf $(OBJ)
