CC = g++-5
CFLAGS = -std=c++11 -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wswitch-default -Wundef -Wno-unused -Wno-unused-parameter
TARGET = main


all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -O2 -o exchange RNG.cpp $(TARGET).cpp  -ltbb -lconfig++ 

debug:
	$(CC) $(CFLAGS) -O0 -g -o exchange RNG.cpp $(TARGET).cpp -ltbb -lconfig++ -pg

clean: 
	$(RM) exchange
