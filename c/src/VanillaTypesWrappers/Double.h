#ifndef _DOUBLE_H_
#define _DOUBLE_H_

struct cDouble
{
    double value;
};
typedef struct cDouble *Double;

Double new_Double(double value);
void del_Double(Double *self);
Double copy_Double(Double self);

#endif // !_DOUBLE_H_