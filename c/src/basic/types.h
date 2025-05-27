#ifndef _TYPES_H_
#define _TYPES_H_

struct cList;
typedef struct cList *List;

struct cString;
typedef struct cString *String;

typedef void (*delfunc)(void **);
typedef void *(*copyfunc)(void *);
typedef String (*strfunc)(void *);

#endif // !_TYPES_H_