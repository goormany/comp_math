CXX = g++
CXXFLAGS = -std=c++17 -O2

all: build run plot

build: main

main: main.cpp
	$(CXX) $(CXXFLAGS) -o main main.cpp

run: main
	./main

plot: run
	python3 ../modules/graphs.py

clean:
	rm -f main nodes.csv interpolation_data.csv interpolation_plot.png