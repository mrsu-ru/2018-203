CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
    main.cpp \
	ivanovii.cpp \
	zhalninrv.cpp \
	scherbakovdv.cpp \
	polyakovda.cpp \
	bagapovar.cpp \
	tarasovams.cpp \
	lab.cpp


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vvm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o vvm
