/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_YY_AMBERMASK_TAB_H_INCLUDED
# define YY_YY_AMBERMASK_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    STRING = 258,
    INUMBER = 259,
    RNUMBER = 260,
    STAR = 261,
    COMMA = 262,
    RANGE = 263,
    RSELECTOR = 264,
    ASELECTOR = 265,
    TSELECTOR = 266,
    RLT = 267,
    RGT = 268,
    ALT = 269,
    AGT = 270,
    NOT = 271,
    AND = 272,
    OR = 273,
    RBRA = 274,
    LBRA = 275,
    ERROR = 276,
    ORIGIN = 277,
    CBOX = 278,
    LIST = 279,
    COM = 280,
    PLANE = 281
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 16 "AmberMask.y" /* yacc.c:1909  */

    SGrammar                gValue;
    SInteger                iValue;
    SReal                   rValue;
    SString                 sValue;
    struct SListItem*       itemValue;
    struct SList*           listValue;
    struct SSelection*      selValue;
    struct SExpression*     exprValue;
    

#line 93 "AmberMask.tab.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_AMBERMASK_TAB_H_INCLUDED  */
