mkdir bin
mkdir bin\shared
mkdir bin\static

mkdir .\bin\static\basic
gcc -c .\src\basic\List\List.c -o .\bin\static\basic\List.o
gcc -c .\src\basic\String\String.c -o .\bin\static\basic\String.o
ar -rcs .\bin\static\libbasic.a .\bin\static\basic\*.o

mkdir bin\static\VanillaTypesWrappers
gcc -c .\src\VanillaTypesWrappers\Double.c -o .\bin\static\VanillaTypesWrappers\Double.o
gcc -c .\src\VanillaTypesWrappers\SizeT.c -o .\bin\static\VanillaTypesWrappers\SizeT.o
ar -rcs .\bin\static\libwrappers.a .\bin\static\VanillaTypesWrappers\*.o

mkdir bin\static\Matrix
gcc -c .\src\Matrix\Matrix.c -o .\bin\static\Matrix\Matrix.o
ar -rcs .\bin\static\libmatrix.a .\bin\static\Matrix\*.o

mkdir bin\static\FLSA_2D
gcc -c .\src\FLSA_2D\FLSA_2D.c -o .\bin\static\FLSA_2D\FLSA_2D.o
ar -rcs .\bin\static\libflsa_2d.a .\bin\static\FLSA_2D\*.o

gcc .\src\main.c -Lbin\static -lflsa_2d -lmatrix -lwrappers -lbasic -o .\bin\static\main.exe