all: tester
	gcc tester.c -o tester number_field.o

number_field.o:
	gcc -c number_field.c number_field.h
