# enable debug under default case
ifeq (x$(debug), x)
  debug=no
endif

source = $(wildcard *.cpp)
ifeq ($(debug), no)
  object = $(patsubst %.cpp, %.g.o, $(source))
  CXXFLAGS = -g -DDEBUG
  lib = -ldeal_II -ltbb -ldl -lpthread $(EXTRA_LIBDIR)
#  LDFLAGS= -ltbb
else
  object = $(patsubst %.cpp, %.o, $(source))
  CXXFLAGS = -g
  lib = -ldeal_II -ltbb -ldl -lpthread $(EXTRA_LIBDIR)
#   LDFLAGS= -ltbb
endif

CXXFLAGS += -I /usr/local/include
LDFLAGS += -L /usr/local/lib

all : main

%.g.o : %.cpp
	$(CXX) -std=c++14 -c -o $@ $< $(CXXFLAGS) $(EXTRA_INCDIR)

%.o : %.cpp
	$(CXX) -std=c++14 -c -o $@ $< $(ARCHS_STANDARD) $(CXXFLAGS) $(EXTRA_INCDIR)

main : $(object)
	$(CXX)  -std=c++14 -o $@ $(object) $(LDFLAGS) $(lib)

clean :
	-rm -rf $(object)
	-rm -rf main
txtclean :
	-rm -rf *.txt

.PHONY : default clean
