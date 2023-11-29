COMP := g++
FLAGS := -std=c++20 -Wall -g

BIN := bin/result.out
SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst src/%.cpp, build/%.o, $(SRC))
INC := -I./src

$(BIN): $(OBJ) 
	$(COMP) $(FLAGS) $(INC) $^ -o $@

build/%.o: src/%.cpp
	@mkdir -p build
	$(COMP) $(FLAGS) -c -o $@ $< $(INC)

clean:
	rm -rf build $(BIN)
	@echo "Cleaned"
