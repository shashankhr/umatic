%{
#include <string.h>
#include "make_structure_reader.h"
void skipcomments(void) ;
extern int endflag;
%}
%option stack
%s IN_STRUCT
%s DEF_STRUCT
%s IN_LINE
%s END_STR
EXPONENT ([Ee][+-]?[0-9]{1,3})
MANTISSA ([0-9]+"."?[0-9]*)|([0-9]*"."?[0-9]+) 
IFDEF "#".*
NUMBER {MANTISSA}{EXPONENT}?
INDEX "["[a-zA-Z0-9_-]+"]"
ARRAY {WORD}{INDEX}+
POINTER "*"+{WSP}*{WORD}
WORD [a-zA-Z][a-zA-Z0-9_-]*
WSP   [ \t]*
OB "\{"
CB "\}"
BEGSTR (typedef{WSP}struct{WSP}?)
ENDSTR ({CB}{WSP}?)
%% 

"exit"  {      
         printf("yylex:Exiting\n");
         exit(0);
        }


<<EOF>> { endflag = 1; return(0);}

{IFDEF}  { printf("yylex:an ifdef %s\n",yytext);
           return(0);
         }
{BEGSTR} {
            BEGIN (DEF_STRUCT);
            printf("yylex:Start struct %s\n",yytext);
            yylval.name = strdup(yytext);
            return(STRSTART);
         }

<IN_STRUCT>{ENDSTR}  {
                        printf("yylex: end of struct\n");
                        BEGIN(END_STR);
                     }

<DEF_STRUCT>{WORD}  {
                      printf("yylex: struct name %s\n",yytext);
                      BEGIN(IN_STRUCT);
                      yylval.name = strdup(yytext);
                      return(STRNAME);
                    }



{POINTER} {
          printf("yylex:a pointer %s\n" ,yytext);
             yylval.name = strdup(yytext);
          if (YYSTATE == IN_LINE){
             return(POINTER);
          }
          printf("yylex: continuing ...\n");
       }

<IN_LINE>{INDEX} { printf("yylex:an index %s\n",yytext);
           yytext[strlen(yytext)-1]=0;
           yylval.index = atoi(yytext+1);
           return(INDEX);
         }

{WORD} {
          printf("yylex:a word %s\n" ,yytext);
             yylval.name = strdup(yytext);
          if (YYSTATE == IN_STRUCT){
             yy_push_state(IN_LINE);
             return(KEY);
          }else if (YYSTATE == IN_LINE){
             return(NAME);
          }else if (YYSTATE == END_STR){
             BEGIN(INITIAL);
             return(STRNAME);
          }else{
             printf("yylex: continuing ...\n");
          }
       }

"/*"   {        
         printf("yylex:Skipping comments\n");
         skipcomments();
       }

{WSP} ;/* do nothing */

"\n" { printf("yylex: end of line\n");
           if(YYSTATE == IN_LINE){
              yy_pop_state();
           }
           return(0);
         }

.  {
    printf("yylex:char: %s\n",yytext);
    yylval.name=strdup(yytext);
    return (OTHER);
   }


%%

void skipcomments(void) {
   char I,O;
   for(
         I=(char)0,O=(char)0;
         ( !((O=='*') && (I=='/')));
         O=I,I=input())
   {
      /*
      if (I=='\n'){
         printf("yylex:End of line comment\n");
         return;
      }
      */
   }
   printf("yylex:Comment delimiter found\n");
   return;
   }

/*
int main(){
  yylex();
  return(0);
}
*/
