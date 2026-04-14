CXX = g++
CXXFLAGS = -Wall -std=c++20

TARGET_MAIN = lb1_main
TARGET_RESEARCH = lb1_research

SRC_MAIN = main.cpp
SRC_RESEARCH = ../modules/research.cpp

all: $(TARGET_MAIN) $(TARGET_RESEARCH)

$(TARGET_MAIN): $(SRC_MAIN)
	$(CXX) $(CXXFLAGS) $(SRC_MAIN) -o $(TARGET_MAIN)

$(TARGET_RESEARCH): $(SRC_RESEARCH)
	$(CXX) $(CXXFLAGS) $(SRC_RESEARCH) -o $(TARGET_RESEARCH)

main: $(TARGET_MAIN)
	./$(TARGET_MAIN)

research: $(TARGET_RESEARCH)
	./$(TARGET_RESEARCH) && python3 ../modules/graphs.py

clean:
	rm -f $(TARGET_MAIN) $(TARGET_RESEARCH) *.csv *.png