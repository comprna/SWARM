EXEC=check_RNA_kit SWARM_preprocess SWARM_preprocess_targets SWARM_preprocess_target_9mers
CC=g++


CFLAGS=-fno-strict-aliasing -DLINUX_ -DOC_NEW_STYLE_INCLUDES -Wno-deprecated -I./htslib-1.19/ -I./slow5lib/include -pthread -D_REENTRANT

LDFLAGS=-L./htslib-1.19/ -L./slow5lib/lib/  -Wl,-R${PWD}/htslib-1.19/ -Wl,-R${PWD}/slow5lib/lib/ -lhts -lslow5

all: $(EXEC)

check_RNA_kit: check_RNA_kit.cpp
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) -g $<

SWARM_preprocess: SWARM_preprocess.cpp
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) -g $<

SWARM_preprocess_targets: SWARM_preprocess_targets.cpp
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) -g $<

SWARM_preprocess_target_9mers: SWARM_preprocess_target_9mers.cpp
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) -g $<

clean: 
	rm -f *.o
	rm -f $(EXEC)
	
.PHONY: all clean
