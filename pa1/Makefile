CPP=g++
CPPFLAGS=-std=c++17 -pedantic -Wvla -Wall -Wshadow -O3
# The build target 
TARGET=pa1

$(TARGET): *.cpp
	$(CPP) $(CPPFLAGS) *.cpp -o $(TARGET) -lm

clean: 
	rm -rf $(TARGET)