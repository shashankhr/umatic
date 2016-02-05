#ifndef SAFEOPEN_H
#define SAFEOPEN_H
extern FILE * safeopen(const char * fname, const char * type);
#ifdef fopen
#undef fopen
#endif
#define fopen safeopen
#endif
