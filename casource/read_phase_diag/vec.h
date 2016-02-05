// vec.h: interface for the vec class.
//
//////////////////////////////////////////////////////////////////////


#ifndef VEC_H
#define VEC_H
#define MAXLEN 1000

typedef double FLOAT;

#include "../machine.h"	
class Mat{
private:
   int my_num;
   int ncol;
   int nrow;
   int ndata;
   int isSquare;
   FLOAT *data;
public:
   Mat();
   Mat(int nc);
   Mat(int nr, int nc);
   Mat(const Mat&);
   virtual ~Mat();
   void init_Mat(int rows,int cols);
   void set_Mat(FLOAT * indata);
   Mat  operator +(Mat  am);
   Mat  operator *(FLOAT sc);
   Mat  operator *(Mat  am);
   const Mat & operator =(const Mat &  am);
   void printMat(FILE * fp);
};


class vec{
private:
	FLOAT x; /* cartesian representation x */
	FLOAT y; /* cartesian representation y */
	FLOAT z;
public:
	vec();
	vec(FLOAT inx,FLOAT iny,FLOAT inz);
	void print(FILE * fp);
	virtual ~vec();
	
	void list (FILE *fp);
	char * strlist();

	void zerovec();
	FLOAT get_x(){ return (x);};
	FLOAT get_y(){ return (y);};
	FLOAT get_z(){ return (z);};
	vec operator +(vec av);
	vec operator -(vec av);
	vec operator * (FLOAT mul);//scalar product
	FLOAT operator *(vec dp);//dot product
	vec operator ^(vec dp);//cross product
	FLOAT operator !(); // magnitude
};
#endif

