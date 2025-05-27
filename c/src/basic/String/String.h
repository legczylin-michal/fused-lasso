#ifndef _STRING_H_
#define _STRING_H_

#include "../types.h"

#include <stdlib.h>
#include <stdbool.h>

String new_String(size_t pLength);
void del_String(String *self);
String copy_String(String self);
String str_String(String self);

size_t String_GetLength(String self);
char *String_GetCharacters(String self);

void String_Append(String self, String pNew);

String String_SubString(String self, int pStart, int pLength);
List String_Split(String pString, String pDelimiter);

String String_FromInt(int value);
String String_FromUnsignedInt(unsigned int value);
String String_FromLongInt(long int value);
String String_FromLongLongInt(long long int value);
String String_FromFloat(float value);
String String_FromDouble(double value);
String String_FromChar(char value);
String String_FromChars(const char *value);
String String_FromBool(bool value);

#define str(X) _Generic((X),               \
    int: String_FromInt,                   \
    unsigned int: String_FromUnsignedInt,  \
    long int: String_FromLongInt,          \
    long long int: String_FromLongLongInt, \
    float: String_FromFloat,               \
    double: String_FromDouble,             \
    char: String_FromChar,                 \
    const char: String_FromChar,           \
    char *: String_FromChars,              \
    const char *: String_FromChars,        \
    bool: String_FromBool)(X)

#endif // !_STRING_H_