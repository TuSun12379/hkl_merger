libsginfo.so: process_hkl.c sg_call.c sgclib.c sgio.c sgfind.c sghkl.c sgsi.c
	gcc -fPIC -shared  process_hkl.c sg_call.c sgclib.c sgio.c sgfind.c sghkl.c sgsi.c -lm -lgsl -lgslcblas -ldl -o libsginfo.so

clean: process_hkl.o sg_call.o sgclib.o sgio.o sgfind.o sghkl.o sgsi.o
	rm -f process_hkl.o sg_call.o sgclib.o sgio.o sgfind.o sghkl.o sgsi.o   