OBJDIR = Object
SRCDIR = Source
INCDIR = Include

HEADERS = $(wildcard $(INCDIR)/*.h)
OBJECTS = $(subst $(INCDIR), $(OBJDIR), $(HEADERS:.h=.o))


CC = g++
CFLAGS = -O2 -std=c++11 `root-config --libs --cflags` \
		 -I$(INCDIR)

LDFLAGS = `root-config --glibs` 

all: $(OBJECTS) analyze_testbeam 
testbeam: $(OBJECTS) analyze_testbeam

analyze_testbeam: $(OBJECTS) analyze_testbeam.cpp 
	@echo "Building executable 'analyze_testbeam'..."
	@$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/Detector.o: $(SRCDIR)/Detector.cpp $(INCDIR)/Detector.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/DetectorSetup.o: $(SRCDIR)/DetectorSetup.cpp $(INCDIR)/DetectorSetup.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@


$(OBJDIR)/TestBeamSetup.o: $(SRCDIR)/TestBeamSetup.cpp $(INCDIR)/TestBeamSetup.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@


$(OBJDIR)/FileReader.o: $(SRCDIR)/FileReader.cpp $(INCDIR)/FileReader.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/TRC_FileReader.o: $(SRCDIR)/TRC_FileReader.cpp $(INCDIR)/TRC_FileReader.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/TRC_channel.o: $(SRCDIR)/TRC_channel.cpp $(INCDIR)/TRC_channel.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@



$(OBJDIR)/TrackInfo.o: $(SRCDIR)/TrackInfo.cpp $(INCDIR)/TrackInfo.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/tracks.o: $(SRCDIR)/tracks.C $(INCDIR)/tracks.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@

$(OBJDIR)/AverageTool.o: $(SRCDIR)/AverageTool.cpp $(INCDIR)/AverageTool.h 
	@echo "Building "$@"..."
	@$(CC) -c $(CFLAGS) $< -o $@


clean: 
	rm $(OBJDIR)/*.o
