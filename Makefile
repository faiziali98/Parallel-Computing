include ./arch.gnu
# OPTIMIZATION = -fast
# OPTIMIZATION = -O3
# DEBUG += -g

app:		cardiacsim2D cardiacsim1D cardiacsimSerial

OBJECTSSerial = cardiacsimSerial.o splot.o cmdLine.o
OBJECTS1D = cardiacsim1D.o splot.o cmdLine.o
OBJECTS2D = cardiacsim2D.o splot.o cmdLine.o

cardiacsimSerial:	        $(OBJECTSSerial) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTSSerial)  $(LDLIBS) -fopenmp
cardiacsim1D:	        $(OBJECTS1D) 
		$(C++LINK) -fopenmp $(LDFLAGS) -o $@ $(OBJECTS1D)  $(LDLIBS)
cardiacsim2D:	        $(OBJECTS2D) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS2D)  $(LDLIBS) -fopenmp

clean:	
	$(RM) *.o cardiacsim1D cardiacsimSerial cardiacsim2D;
	$(RM) core;
