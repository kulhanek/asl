/* =============================================================================
AMBER Mask Lexical Analyzer
      Copyright (c) 2008 Petr Kulhanek, kulhanek@enzim.hu
============================================================================= */

%{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "AmberMaskParser.hpp"
#include "AmberMask.tab.h"
%}

%option noyywrap
newline       [\n]
whitesp       [ \t]+
string        [0-9]*[a-zA-Z]+[a-zA-Z0-9$^+_*=\-']*
rnumber       [0-9]+[.]?[0-9]*
inumber       [0-9]+

/* rules ==================================================================== */
%%

"origin"    { yylval.gValue.Position = LexPosition; LexPosition += 6; return ORIGIN; }
"cbox"      { yylval.gValue.Position = LexPosition; LexPosition += 4; return CBOX; }
"list"      { yylval.gValue.Position = LexPosition; LexPosition += 4; return LIST; }
"com"       { yylval.gValue.Position = LexPosition; LexPosition += 3; return COM; }
"plane"     { yylval.gValue.Position = LexPosition; LexPosition += 5; return PLANE; }

"*"       { yylval.gValue.Position = LexPosition; LexPosition++; return STAR; }
","       { yylval.gValue.Position = LexPosition; LexPosition++; return COMMA; }
"-"       { yylval.gValue.Position = LexPosition; LexPosition++; return RANGE; }

":"       { yylval.gValue.Position = LexPosition; LexPosition++; return RSELECTOR; }
"@"       { yylval.gValue.Position = LexPosition; LexPosition++; return ASELECTOR; }
"@%"      { yylval.gValue.Position = LexPosition; LexPosition += 2; return TSELECTOR; }

"<:"      { yylval.gValue.Position = LexPosition; LexPosition += 2; return RLT; }
">:"      { yylval.gValue.Position = LexPosition; LexPosition += 2; return RGT; }
"<@"      { yylval.gValue.Position = LexPosition; LexPosition += 2; return ALT; }
">@"      { yylval.gValue.Position = LexPosition; LexPosition += 2; return AGT; }

"&"       { yylval.gValue.Position = LexPosition; LexPosition++; return AND; }
"|"       { yylval.gValue.Position = LexPosition; LexPosition++; return OR; }
"!"       { yylval.gValue.Position = LexPosition; LexPosition++; return NOT; }

"("       { yylval.gValue.Position = LexPosition; LexPosition++; return LBRA; }
")"       { yylval.gValue.Position = LexPosition; LexPosition++; return RBRA; }

{string}  {
    yylval.sValue.Position = LexPosition;
    if( yyleng > 4 ){
        pperror("name is too long (max 4 characters)",LexPosition);
        }
    memset(yylval.sValue.String,' ',4);
    strncpy(yylval.sValue.String,yytext,yyleng);
    LexPosition += yyleng;
    return STRING;
    }

{inumber} {
    yylval.iValue.Position = LexPosition;
    if( sscanf(yytext, "%d", &yylval.iValue.Number) != 1 ){
        pperror("unable to convert string to integer",LexPosition);
        }
    LexPosition += yyleng;
    return INUMBER;
    }

{rnumber} {
    yylval.rValue.Position = LexPosition;
    if( sscanf(yytext, "%lf", &yylval.rValue.Number) != 1 ){
        pperror("unable to convert string to real number",LexPosition);
        }
    return RNUMBER;
    }

{whitesp} { LexPosition += yyleng;}

{newline} {
            pperror("new line character is not allowed in mask specification",LexPosition);
            return ERROR;
            }
