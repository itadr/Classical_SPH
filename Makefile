TARGET= main
CXX= g++
CXXFLAGS= -std=c++17
SUFFIX= .cpp
SRCS= $(wildcard *$(SUFFIX))
OBJS= $(SRCS:$(SUFFIX)=.o)

$(TARGET):$(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) 
%.o:%.cpp
	$(CXX) -c $^ $(CXXFLAGS)
clean:
	rm -f *.o main *.csv 