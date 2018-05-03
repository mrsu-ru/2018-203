CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
    main.cpp \
	ivanovii.cpp \
	ivanovdd.cpp \
	zhalninrv.cpp \
	scherbakovdv.cpp \
	polyakovda.cpp \
	grishaevaoov.cpp \
	bagapovar.cpp \
	tarasovams.cpp \
	syusinaev.cpp\
	fedyanovaam.cpp \
	serguninaes.cpp \
	borisovrs.cpp \
	kuznetsovais.cpp\
	nefedovms.cpp \
	bulychevaoa.cpp \
	salinaa.cpp \
	seninvs.cpp \
	mescheryakovam.cpp \
	itaevde.cpp \
	biryukovaes.cpp \
	velmiskinaav.cpp \
	kozlovdn.cpp\
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
