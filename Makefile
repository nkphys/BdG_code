OBJS = functions_real.o main.o

EXENAME  = BdG
### ------ Personal PC compilation ------------
OPTFLAG = -O3
CC = g++ $(OPTFLAG) -std=c++11
CFLAGS = -c
#-DUSE_COMPLEX

LDFLAGS  = -llapack -lblas

### --- turn on for debugging -----------
CPPFLAGS += -Isrc
CPPFLAGS += -O3 #### Reduce compilation time and make debugging produce the expected results.
STRIP_COMMAND = true #### Keeps lines in the executable for debugging

$(EXENAME): clean $(OBJS) 
	$(CC) $(OBJS) -o $(EXENAME) $(LDFLAGS)
	$(STRIP_COMMAND) $(EXENAME)

all:	$(EXENAME)
	 
clean:
	rm -f $(EXENAME) *.o


functions_real.o : functions_real.cpp
	$(CC) $(CFLAGS) functions_real.cpp $(LDFLAGS)

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp $(LDFLAGS)


######## End of Makefile ########
