CXX = g++
CXXFLAGS = -O2 -Wall
LDFLAGS = -llapack -lblas

GEN_SRC = ../modules/generator.cpp
SOLVER_SRC = main.cpp

GEN_EXE = generator
SOLVER_EXE = solver

all: $(GEN_EXE) $(SOLVER_EXE)

$(GEN_EXE): $(GEN_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(SOLVER_EXE): $(SOLVER_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

run-gen: $(GEN_EXE)
	./$(GEN_EXE)

run: $(SOLVER_EXE)
	./$(SOLVER_EXE)

clean:
	rm -f $(GEN_EXE) $(SOLVER_EXE) input.txt solution.txt

.PHONY: all run-gen run clean