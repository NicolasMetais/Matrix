NAME = Matrix
CC = c++
CFLAGS = -Wall -Wextra -Werror -std=c++17 -g3

SRC_DIR = ex
EXS = $(shell seq -w 0 14)

EX_TARGETS = $(addprefix $(SRC_DIR), $(addsuffix /$(NAME), $(EXS)))

all: $(EX_TARGETS)

%/Matrix:
	@echo "Compiling $*..."
	$(CC) $(CFLAGS) $*/main.cpp -o  $@

clean:
	rm -rf $(EX_TARGETS)

re: clean all

.PHONY: all clean re
