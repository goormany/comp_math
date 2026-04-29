CXX = g++
CXXFLAGS = -O2 -Wall
LDFLAGS = -llapack -lblas

SOLVER_SRC = main.cpp

SOLVER_EXE = solver

all: $(SOLVER_EXE)

$(SOLVER_EXE): $(SOLVER_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

run: $(SOLVER_EXE)
	./$(SOLVER_EXE)

clean:
	rm -f $(SOLVER_EXE)

.PHONY: all run clean