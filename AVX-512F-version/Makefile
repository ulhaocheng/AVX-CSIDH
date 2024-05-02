# compiler 
CC = gcc-9

# include folder 
INC_DIR += -I ./src/inc

# the required files for batched high-throughput (ht) implementations
HT_FILES =  ./src/lib/rdtsc64.S \
            ./src/lib/fp512.S \
            ./src/lib/gfparith.c \
            ./src/lib/tedcurve.c \
            ./src/lib/rng.c \
            ./src/lib/utils.c \
            ./src/lib/main_ht.c 

# the specific action file and header file
HT_ACTION_FILE = ./src/action/$(shell echo $(STYLE))/action_$(shell echo $(METHOD) | tr A-Z a-z).c
ACTION_INC = -I ./src/inc/action_inc/$(shell echo $(STYLE))

# flags
FLAGS = -O2 -mavx512f -fomit-frame-pointer -m64 -mbmi2 -march=native -fwrapv -mtune=native

help:
	@echo "\nCompile the batched high-throughput (ht) implementation:"
	@echo "make ht_csidh STYLE=[OAYT/DUMMY_FREE] METHOD=[EXTRA_DUMMY/EXTRA_INFINITY/COMBINED]"
	@echo "\nClean the exectuable files:"
	@echo "make clean"

ht_csidh: 
	$(CC) $(INC_DIR) $(ACTION_INC) $(HT_FILES) $(HT_ACTION_FILE) -o ht_csidh $(FLAGS)

clean:
	rm -f ht_csidh 