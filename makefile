
OBJS = SimpleTectonics.cpp

TINYLINK = -lTinyEngine -lX11 -lpthread -lSDL2 -lnoise -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lGL -lGLEW -lboost_serialization -lboost_system -lboost_filesystem

CC = g++ -std=c++17
COMPILER_FLAGS = -O
LINKER_FLAGS = -I/usr/local/include -L/usr/local/lib
OBJ_NAME = tectonics

all: $(OBJS)
			$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) $(TINYLINK) -o $(OBJ_NAME)
