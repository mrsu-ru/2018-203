CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
<<<<<<< HEAD
        main.cpp \
		ivanovdd.cpp \
				ivanovii.cpp \
				lab.cpp

=======
    main.cpp \
	ivanovii.cpp \
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
	mescheryakovam.cpp \
	itaevde.cpp \
	biryukovaes.cpp \
	lab.cpp
>>>>>>> 96a988efb4ec7cf7cd7777aef93ef755bc38c6e3

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vvm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o vvm
