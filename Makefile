TARGET = test.exe
OBJS = main.o 
CC = g++
CFLAGS = -c -Wall  -g -std=c++11
LFLAGS = -Wall  -g -std=c++11
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFALGS = -Wall -O3 -DNDEBUG -std=c++11

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp xy.h vfield.h \
        system.h density.h orientation.h flux.h
	$(CC) $(CFLAGS) main.cpp


.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: cleandata
cleandata:
	rm -f *.dat
	rm -f results/*.dat
  

