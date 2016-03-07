CC = g++-4.8
CFLAGS = -std=c++11 -pthread -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused
TARGET = main


all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -O2 -o exchange RNG.cpp $(TARGET).cpp -lconfig++ -lgflags 
#if I need to use boost libraries, I need to fix the link with:
#-I/usr/local/gcc/include -L/usr/local/gcc/lib -lboost_regex
# install_name_tool -change libboost_regex.dylib /usr/local/gcc/lib/libboost_regex.dylib exchange

debug:
	$(CC) $(CFLAGS) -O0 -g -o exchange RNG.cpp $(TARGET).cpp -lconfig++ -lgflags 

clean: 
	$(RM) exchange