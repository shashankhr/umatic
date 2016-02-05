#!/usr/local/bin/gnumake
vpath %.c  ..
#vpath %.h  ..
#vpath %.o  ..

.y.c:
	$(YACC) -d -o $@ $<

CXX=g++

YACC = bison
LEX = flex
.PHONY: res
.PHONY: fsres
.PHONY: allres
.PHONY: struct
.PHONY: mould

#define the DIRT to be removed by make clean
DIRT=\
P* \
W* \
_*\
make_structure_reader.[ch]\
make_structure_output.[ch]\
make_structure_reader\
make_structure_output\
i_res\
f_res_single\
f_res_double\
fsres_single\
fsres_double\
Struct.in\
safeopen.o\
makemould.o\
mould.inp\
mould.png\
makemould\
mdeb


default: res

# the dependencies
val: write_block_values.c write_ctrl_values.c
struct: make_structure_output make_structure_reader
res: f_res_single f_res_double i_res
fsres: fsres_single fsres_double
allres: res fsres

f_res_double f_res_single: f_res.c 
fsres_double fsres_single: fsres.c 
i_res : i_res.c 

fsres_double fsres_single f_res_double f_res_single i_res : debin.c safeopen.o

f_res_double fsres_double i_res:
	$(CC) -DCA_DOUBLE -I .. -I/usr/local/include $^  -o $@  -lm -L/usr/local/lib -lgd  -lpng -lz   
f_res_single fsres_single:
	$(CC) -I .. -I/usr/local/include $^  -o $@  -lm -L/usr/local/lib -lgd -lpng -lz  

write_block_values.c: readstruct blocks.h
	readstruct < blocks.h > lexout.txt;mv P_out.c $@;mv P_out.h $*.h

write_ctrl_values.c: readstruct read_ctrl.h
	readstruct < read_ctrl.h > lexout.txt;mv P_out.c $@;mv P_out.h $*.h

write_matprop_values.c: readstruct matprops.h
	readstruct < matprops.h > lexout.txt;mv P_out.c $@;mv P_out.h $*.h

readlexer.c: make_structure_reader.c

make_structure_output:	make_structure_output.c readlexer.c make_structure_reader.h
	cc -g  make_structure_output.c readlexer.c -o make_structure_output -ll

make_structure_reader:	make_structure_reader.c readlexer.c make_structure_reader.h
	cc -g  make_structure_reader.c readlexer.c -o make_structure_reader -ll
clean:
	rm -vf $(DIRT)

#the makemould program
mould: makemould
	./makemould
makemould.o:	makemould.cxx
	$(CXX) -c $<
makemould:  makemould.o
	$(CXX) -o $@ $<  -L/usr/local/lib -lgd -lm -lz
mdeb:  makemould.cxx
	$(CXX) -ggdb -o $@ $<  -lefence -L/usr/local/lib -lgd -lm -lpng -lz
