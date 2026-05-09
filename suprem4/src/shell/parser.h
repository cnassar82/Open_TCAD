/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    PARAMETER = 258,               /* PARAMETER  */
    COMMAND = 259,                 /* COMMAND  */
    NAME = 260,                    /* NAME  */
    LIST = 261,                    /* LIST  */
    SOURCE = 262,                  /* SOURCE  */
    EOL = 263,                     /* EOL  */
    ENDFILE = 264,                 /* ENDFILE  */
    QUIT = 265,                    /* QUIT  */
    BACK = 266,                    /* BACK  */
    REDIRECT = 267,                /* REDIRECT  */
    BANG = 268,                    /* BANG  */
    HELP = 269,                    /* HELP  */
    BG_GRP = 270,                  /* BG_GRP  */
    END_GRP = 271,                 /* END_GRP  */
    DELIMIT = 272,                 /* DELIMIT  */
    DEFINE = 273,                  /* DEFINE  */
    COMMENT = 274,                 /* COMMENT  */
    UNDEF = 275,                   /* UNDEF  */
    SET = 276,                     /* SET  */
    UNSET = 277,                   /* UNSET  */
    FOR = 278,                     /* FOR  */
    END = 279,                     /* END  */
    CHDIR = 280                    /* CHDIR  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif
/* Token kinds.  */
#define YYEMPTY -2
#define YYEOF 0
#define YYerror 256
#define YYUNDEF 257
#define PARAMETER 258
#define COMMAND 259
#define NAME 260
#define LIST 261
#define SOURCE 262
#define EOL 263
#define ENDFILE 264
#define QUIT 265
#define BACK 266
#define REDIRECT 267
#define BANG 268
#define HELP 269
#define BG_GRP 270
#define END_GRP 271
#define DELIMIT 272
#define DEFINE 273
#define COMMENT 274
#define UNDEF 275
#define SET 276
#define UNSET 277
#define FOR 278
#define END 279
#define CHDIR 280

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 56 "shell/parser.y"

    int ival;
    char *sval;
    

#line 123 "y.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
