#include "String.h"

#include <stdio.h>

#include "../List/List.h"

struct cString
{
    size_t _length;
    char *_characters;
};

String new_String(size_t pLength)
{
    String result = (String)malloc(sizeof(struct cString));

    result->_length = pLength;
    result->_characters = (char *)malloc(sizeof(char) * (result->_length + 1));
    result->_characters[result->_length] = '\0';

    return result;
}

void del_String(String *self)
{
    if (*self == NULL)
        return;

    free((*self)->_characters);
    free(*self);

    *self = NULL;

    return;
}

String copy_String(String self)
{
    String result = new_String(self->_length);

    for (size_t i = 0; i < result->_length; i++)
        result->_characters[i] = self->_characters[i];

    return result;
}

String str_String(String self)
{
    return copy_String(self);
}

size_t String_GetLength(String self)
{
    return self->_length;
}

char *String_GetCharacters(String self)
{
    return self->_characters;
}

void String_Append(String self, String pNew)
{
    self->_characters = realloc(self->_characters, sizeof(char) * (self->_length + pNew->_length + 1));
    for (size_t i = 0; i < pNew->_length; i++)
        self->_characters[self->_length + i] = pNew->_characters[i];
    self->_length = self->_length + pNew->_length;
    self->_characters[self->_length] = '\0';
}

String String_SubString(String self, int pStart, int pLength)
{
    size_t _start = pStart < 0 ? self->_length - ((-pStart) % self->_length) : pStart % self->_length;
    int direction = pLength < 0 ? -1 : 1;
    size_t _length = pLength < 0 ? -pLength : pLength;
    if (direction == 1 && _start + _length > self->_length)
        _length = self->_length - _start;
    else if (direction == -1 && (int)_start - (int)_length < 0)
        _length = _start + 1;

    String result = new_String(_length);

    for (size_t i = 0; i < _length; i++)
        result->_characters[i] = self->_characters[_start + direction * i];

    return result;
}

List String_Split(String target, String delimiter)
{
    List result = new_List((delfunc)del_String, (copyfunc)copy_String, (strfunc)str_String);

    int from = 0;
    for (size_t i = 0; i < target->_length - delimiter->_length + 1; i++)
    {
        bool delimiterFound = true;
        for (size_t j = 0; j < delimiter->_length; j++)
        {
            if (target->_characters[i + j] == delimiter->_characters[j])
                continue;

            delimiterFound = false;
            break;
        }

        if (delimiterFound)
        {
            int to = i - 1;
            List_Append(result, String_SubString(target, from, to - from + 1));
            from = i + delimiter->_length;
        }
    }
    List_Append(result, String_SubString(target, from, target->_length - from));

    return result;
}

String String_FromInt(int value)
{
    int bufferSize = snprintf(NULL, 0, "%d", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%d", value);

    return result;
}

String String_FromUnsignedInt(unsigned int value)
{
    int bufferSize = snprintf(NULL, 0, "%d", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%d", value);

    return result;
}

String String_FromLongInt(long int value)
{
    int bufferSize = snprintf(NULL, 0, "%ld", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%ld", value);

    return result;
}

String String_FromLongLongInt(long long int value)
{
    int bufferSize = snprintf(NULL, 0, "%lld", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%lld", value);

    return result;
}

String String_FromFloat(float value)
{
    int bufferSize = snprintf(NULL, 0, "%f", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%f", value);

    return result;
}

String String_FromDouble(double value)
{
    int bufferSize = snprintf(NULL, 0, "%lf", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%lf", value);

    return result;
}

String String_FromChar(char value)
{
    int bufferSize = snprintf(NULL, 0, "%c", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%c", value);

    return result;
}

String String_FromChars(const char *value)
{
    int bufferSize = snprintf(NULL, 0, "%s", value);

    String result = new_String(bufferSize);

    snprintf(result->_characters, bufferSize + 1, "%s", value);

    return result;
}

String String_FromBool(bool value)
{
    return value ? String_FromChars("true") : String_FromChars("false");
}