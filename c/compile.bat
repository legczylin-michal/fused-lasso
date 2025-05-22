mkdir bin
mkdir bin\shared
mkdir bin\static

mkdir bin\static\matrix
gcc -c .\src\matrix\Matrix.c -o .\bin\static\matrix\Matrix.o
ar -rcs .\bin\static\libmatrix.a .\bin\static\matrix\*.o

mkdir bin\static\flsa_2d
gcc -c .\src\flsa_2d\FLSA_2D.c -o .\bin\static\flsa_2d\FLSA_2D.o
ar -rcs .\bin\static\libflsa_2d.a .\bin\static\flsa_2d\*.o

gcc .\src\main.c -Lbin\static -lmatrix  -lflsa_2d -o .\bin\static\main.exe