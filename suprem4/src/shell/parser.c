/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "shell/parser.y"


/*************************************************************************
 *									 *
 *   Original : MEL         Stanford University        Sept, 1984	 *
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *     Copyright c 1991 The board of regents of the University of 	 *
 *                      Florida.  All rights reserved.			 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of the University of 	 *
 *     Florida.								 *
 *									 *
 *									 *
 *************************************************************************/
/*   parser.y                Version 5.1     */
/*   Last Modification : 7/3/91 08:28:19 */

/************************************************************************
 *									*
 *	yyparse - this yacc file specifies the grammer used in the 	*
 *  parser.  The return value is -1 for end of input, 2 for continue    *
 *  with a prompt, and 0 for continue with no prompt.			*
 *									*
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "sysdep.h"


static char *tp[10];
static int te[10];

#ifdef MEMDEB
#define malloc MALLOC
#endif

#include "shell.h"

#define YYRETURN(A)        { return(A); }


#line 122 "y.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
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
    

#line 231 "y.tab.c"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_PARAMETER = 3,                  /* PARAMETER  */
  YYSYMBOL_COMMAND = 4,                    /* COMMAND  */
  YYSYMBOL_NAME = 5,                       /* NAME  */
  YYSYMBOL_LIST = 6,                       /* LIST  */
  YYSYMBOL_SOURCE = 7,                     /* SOURCE  */
  YYSYMBOL_EOL = 8,                        /* EOL  */
  YYSYMBOL_ENDFILE = 9,                    /* ENDFILE  */
  YYSYMBOL_QUIT = 10,                      /* QUIT  */
  YYSYMBOL_BACK = 11,                      /* BACK  */
  YYSYMBOL_REDIRECT = 12,                  /* REDIRECT  */
  YYSYMBOL_BANG = 13,                      /* BANG  */
  YYSYMBOL_HELP = 14,                      /* HELP  */
  YYSYMBOL_BG_GRP = 15,                    /* BG_GRP  */
  YYSYMBOL_END_GRP = 16,                   /* END_GRP  */
  YYSYMBOL_DELIMIT = 17,                   /* DELIMIT  */
  YYSYMBOL_DEFINE = 18,                    /* DEFINE  */
  YYSYMBOL_COMMENT = 19,                   /* COMMENT  */
  YYSYMBOL_UNDEF = 20,                     /* UNDEF  */
  YYSYMBOL_SET = 21,                       /* SET  */
  YYSYMBOL_UNSET = 22,                     /* UNSET  */
  YYSYMBOL_FOR = 23,                       /* FOR  */
  YYSYMBOL_END = 24,                       /* END  */
  YYSYMBOL_CHDIR = 25,                     /* CHDIR  */
  YYSYMBOL_YYACCEPT = 26,                  /* $accept  */
  YYSYMBOL_line = 27,                      /* line  */
  YYSYMBOL_commands = 28,                  /* commands  */
  YYSYMBOL_endline = 29,                   /* endline  */
  YYSYMBOL_rdir = 30,                      /* rdir  */
  YYSYMBOL_parameter = 31,                 /* parameter  */
  YYSYMBOL_act1 = 32,                      /* act1  */
  YYSYMBOL_act2 = 33                       /* act2  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int8 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  37
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   132

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  26
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  8
/* YYNRULES -- Number of rules.  */
#define YYNRULES  30
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  68

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   280


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,    73,    73,    98,   117,   145,   163,   183,   228,   287,
     319,   327,   341,   355,   382,   389,   396,   403,   415,   416,
     421,   423,   425,   427,   432,   437,   442,   448,   453,   462,
     486
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "PARAMETER", "COMMAND",
  "NAME", "LIST", "SOURCE", "EOL", "ENDFILE", "QUIT", "BACK", "REDIRECT",
  "BANG", "HELP", "BG_GRP", "END_GRP", "DELIMIT", "DEFINE", "COMMENT",
  "UNDEF", "SET", "UNSET", "FOR", "END", "CHDIR", "$accept", "line",
  "commands", "endline", "rdir", "parameter", "act1", "act2", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-39)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-1)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int8 yypact[] =
{
     107,    -6,     0,     0,   -39,   -39,   -39,     6,     0,     0,
     -39,   -39,     0,     0,     0,     0,     0,     3,     0,    15,
     -39,   -39,   -39,    -8,    -8,   -39,   -39,    -8,    -8,   107,
      -6,    -6,    -6,    -6,    -6,    18,    -6,   -39,    22,    -6,
      -6,    -6,    -6,   -39,    59,   -39,   -39,   -39,   -39,   -39,
     -39,   -39,   -39,   -39,   -39,   -39,   -39,   -39,   -39,   107,
      -8,    82,    -6,   -39,   -39,    -8,    -6,   -39
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       0,     0,    28,    28,    20,    15,    16,     0,    28,    28,
      29,    22,    28,    28,    28,    28,    28,     0,    28,     0,
      14,    17,    27,    26,    26,    21,    23,    26,    26,     0,
       0,     0,     0,     0,     0,     0,     0,     1,    24,     0,
       0,     0,     0,    19,     0,     6,    10,     5,    11,    12,
      29,    13,    25,     2,     3,     7,     4,    30,    18,     0,
      26,     0,     0,    30,     9,    26,     0,     8
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -39,     1,   -38,    -1,   -18,     4,   -24,   -36
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,    43,    44,    20,    39,    23,    29,    60
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      21,    19,     4,    22,    38,     7,    40,    24,    35,    41,
      42,    11,    27,    28,    25,    37,    30,    31,    32,    33,
      34,    61,    36,    26,    50,    52,    59,    65,     0,    45,
      46,    47,    48,    49,     0,    51,     0,     0,    53,    54,
      55,    56,    62,     0,     0,    58,     0,    66,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       1,    64,    58,     2,     0,    67,     3,     4,     5,     6,
       7,     0,     8,     9,    10,    57,    11,    12,    13,    14,
      15,    16,    17,     1,    18,     0,     2,     0,     0,     3,
       4,     5,     6,     7,     0,     8,     9,    10,     0,    11,
      12,    13,    14,    15,    16,    17,    63,    18,     1,     0,
       0,     2,     0,     0,     3,     4,     5,     6,     7,     0,
       8,     9,    10,     0,    11,    12,    13,    14,    15,    16,
      17,     0,    18
};

static const yytype_int8 yycheck[] =
{
       1,     0,     8,     3,    12,    11,    24,     3,     5,    27,
      28,    17,     8,     9,     8,     0,    12,    13,    14,    15,
      16,    59,    18,    17,     6,     3,    50,    63,    -1,    30,
      31,    32,    33,    34,    -1,    36,    -1,    -1,    39,    40,
      41,    42,    60,    -1,    -1,    44,    -1,    65,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       1,    62,    61,     4,    -1,    66,     7,     8,     9,    10,
      11,    -1,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,     1,    25,    -1,     4,    -1,    -1,     7,
       8,     9,    10,    11,    -1,    13,    14,    15,    -1,    17,
      18,    19,    20,    21,    22,    23,    24,    25,     1,    -1,
      -1,     4,    -1,    -1,     7,     8,     9,    10,    11,    -1,
      13,    14,    15,    -1,    17,    18,    19,    20,    21,    22,
      23,    -1,    25
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,     1,     4,     7,     8,     9,    10,    11,    13,    14,
      15,    17,    18,    19,    20,    21,    22,    23,    25,    27,
      29,    29,     3,    31,    31,     8,    17,    31,    31,    32,
      31,    31,    31,    31,    31,     5,    31,     0,    12,    30,
      30,    30,    30,    27,    28,    29,    29,    29,    29,    29,
       6,    29,     3,    29,    29,    29,    29,    16,    27,    32,
      33,    28,    30,    24,    29,    33,    30,    29
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    26,    27,    27,    27,    27,    27,    27,    27,    27,
      27,    27,    27,    27,    27,    27,    27,    27,    28,    28,
      29,    29,    29,    29,    30,    30,    30,    31,    31,    32,
      33
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     4,     4,     4,     3,     3,     4,     9,     7,
       3,     3,     3,     3,     1,     1,     1,     2,     2,     1,
       1,     2,     1,     2,     1,     2,     0,     1,     0,     0,
       0
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* line: COMMAND parameter rdir endline  */
#line 74 "shell/parser.y"
        {
	    char *cmd = (yyvsp[-3].sval), *par = (yyvsp[-2].sval), *red = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);

	    do_echo();

	    /*do no take the action if storing*/
	    if (depth == -1) {
		/*call the command setup routine*/
		if ((par == NULL) || (interactive))
		    do_command(cmd, par, TRUE, red, (back & BACKGROUND) );
		else
		    do_command(cmd, par, FALSE, red, (back & BACKGROUND) );
		}

	    /*free the used malloc space*/
	    free(cmd);
	    if (par != NULL) free(par);
	    if (red != NULL) free(red);

	    if (depth == -1) {
		YYRETURN( back & PROMPT );
	    }
	}
#line 1311 "y.tab.c"
    break;

  case 3: /* line: SOURCE parameter rdir endline  */
#line 99 "shell/parser.y"
        {
	    char *par = (yyvsp[-2].sval), *red = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);

	    do_echo();

	    /*do not take the action if storing*/
	    if (depth == -1) {
		do_source(par, red, ( back & BACKGROUND), TRUE );

	        if ( (back & PROMPT) && (back & BACKGROUND) ) {
		    YYRETURN( PROMPT );
		}
		else {
		    YYRETURN( ~ PROMPT );
		}
	    }
	}
#line 1334 "y.tab.c"
    break;

  case 4: /* line: HELP parameter rdir endline  */
#line 118 "shell/parser.y"
        {
	    char *par = (yyvsp[-2].sval), *red = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);

	    FILE tout, *t1;
	    do_echo();

	    /*do the redirection if necessary*/
	    tout = *stdout;
	    if (red != NULL)  {
		t1 = file_parse(red, "a");
		*stdout = *t1;
	    }

	    help( par );

	    if (red != NULL) {
		fclose(stdout);
		*stdout = tout;
	    }

	    /*free the space*/
	    if (par != NULL) free(par);
	    if (red != NULL) free(red);

	    if (depth == -1) YYRETURN( back & PROMPT );
	}
#line 1366 "y.tab.c"
    break;

  case 5: /* line: UNDEF parameter endline  */
#line 146 "shell/parser.y"
        {
	    char *par = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);

	    do_echo();

	    /*do no take the action if storing*/
	    if (depth == -1) {
		/*do the undefintion*/
		if (par != NULL)
		    if ( undef_macro(par, &macro) == -1)
			fprintf(stderr, "can not find macro %s\n", par);
		}
	    if (par != NULL) free(par);
	    if (depth == -1)
		YYRETURN( back & PROMPT );
	}
#line 1388 "y.tab.c"
    break;

  case 6: /* line: DEFINE parameter endline  */
#line 164 "shell/parser.y"
        {
	    char *par = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);

	    do_echo();

	    /*do no take the action if storing*/
	    if (depth == -1) {
		/*handle the definitions and results*/
		if (par == NULL)
		    dump_macro(macro);
		else
		    define_macro(par, &macro);
		}

	    if (par != NULL) free(par);
	    if (depth == -1)
		YYRETURN( back & PROMPT );
	}
#line 1412 "y.tab.c"
    break;

  case 7: /* line: BANG parameter rdir endline  */
#line 184 "shell/parser.y"
        {
	    char *par = (yyvsp[-2].sval), *red = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);
	    char str[160];
	    int status, pid, w;

	    do_echo();

	    str[0] = '\0';

	    /*do no take the action if storing*/
	    if (depth == -1) {
		/*if no command give a shell*/
		if (par == NULL)
		    strcpy(str, "csh");
		else {
		    strcpy(str, par);
		    free(par);
		    }

		/*if a redirect add that on*/
		if (red != NULL) {
		    strcat(str, red);
		    free(red);
		    }

		/*add background on*/
		if (back & BACKGROUND)
		    strcat(str, " & ");

		/*exec a shell to handle the request*/
		if ((pid = vfork()) == 0) {
		    if (execl("/bin/sh", "sh", "-c", str, 0) == -1)
			printf("error number %d\n", 127);
		    _exit(127);
		}

		/*wait for the process to finish*/
		while ((w = wait(&status)) != pid && w != -1);

	    }
	    if (depth == -1)
		YYRETURN( PROMPT );
	}
#line 1461 "y.tab.c"
    break;

  case 8: /* line: FOR NAME LIST act1 commands END act2 rdir endline  */
#line 229 "shell/parser.y"
        {
	    char *body;
	    char *list;
	    char *name;
	    char *rdir;
	    char *value;
	    char str[160];
	    int endline;
	    int i;
	    extern char *list_parse();

	    echo_bufptr = 0;

	    prompt = tp[depth];
	    echoall = te[depth];
	    depth--;
	    if ( depth == -1 ) {
		/*we need to save the names, because parsing will
		  destroy the current set*/
		body = (yyvsp[-2].sval);
		list = (yyvsp[-6].sval);
		name = (yyvsp[-7].sval);
		rdir = (yyvsp[-1].sval);
		endline = (yyvsp[0].ival);

		/*strip the "end" off*/
		for(i = strlen(body)-1; body[i] != 'e'; i--);
		body[i] = '\0';

		for( value = list_parse( &list, TRUE ); value != NULL;
		     value = list_parse( &list, FALSE) ) {
		    /*basic plan is to use macro processor to define strings*/
		    strcpy(str, name);
		    strcat(str, " ");
		    strcat(str, value);
		    free(value);
		    define_macro(str, &macro);

		    do_string( body, rdir, ( endline & BACKGROUND ) );
		}


		/*remove the entry from the macro table*/
		(void)undef_macro(name, &macro);

		if (body != NULL) free(body);
		if (list != NULL) free(list);
		if (name != NULL) free(name);
		if (rdir != NULL) free(rdir);

		if ( (endline & PROMPT) || (endline & BACKGROUND) ) {
		    YYRETURN( PROMPT );
		}
		else {
		    YYRETURN( ~ PROMPT );
		}
	    }
	}
#line 1524 "y.tab.c"
    break;

  case 9: /* line: BG_GRP act1 commands END_GRP act2 rdir endline  */
#line 288 "shell/parser.y"
        {
	    char *body = (yyvsp[-2].sval), *red = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);
	    char *s;
	    int i;

	    echo_bufptr = 0;

	    prompt = tp[depth];
	    echoall = te[depth];
	    depth--;
	    /*do no take the action if storing*/
	    if (depth == -1)  {
		s = body;

		/*strip off the trailing }*/
		for(i = strlen(s)-1; s[i] != '}'; i--);
		s[i] = ';';

		do_string(body, red, ( back & BACKGROUND) );

		if (body != NULL) free(body);

		if ( (back & PROMPT) || (back & BACKGROUND) ) {
		    YYRETURN( PROMPT );
		}
		else {
		    YYRETURN( ~ PROMPT );
		}
	    }
	}
#line 1560 "y.tab.c"
    break;

  case 10: /* line: COMMENT parameter endline  */
#line 320 "shell/parser.y"
        {
	    int back = (yyvsp[0].ival);
	    do_echo();

	    if (depth == -1)
		YYRETURN( back & PROMPT );
	}
#line 1572 "y.tab.c"
    break;

  case 11: /* line: SET parameter endline  */
#line 328 "shell/parser.y"
        {
	    char *par = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);
	    do_echo();

	    /*do no take the action if storing*/
	    if (depth == -1)
		do_set(par, TRUE);
	    /*free up the wasted space*/
	    if (par != NULL) free(par);
	    if (depth == -1)
		YYRETURN( back & PROMPT );
	}
#line 1590 "y.tab.c"
    break;

  case 12: /* line: UNSET parameter endline  */
#line 342 "shell/parser.y"
        {
	    char *par = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);
	    do_echo();

	    /*do no take the action if storing*/
	    if (depth == -1)
		do_set(par, FALSE);
	    /*free up the wasted space*/
	    if (par != NULL) free(par);
	    if (depth == -1)
		YYRETURN( back & PROMPT );
	}
#line 1608 "y.tab.c"
    break;

  case 13: /* line: CHDIR parameter endline  */
#line 356 "shell/parser.y"
        {
	    char *par = (yyvsp[-1].sval);
	    int back = (yyvsp[0].ival);
	    char *s;

	    do_echo();

	    /*do not take the action if storing*/
	    if (depth == -1) {
		s = par;
		if ( s == NULL ) {
		    if (chdir( getenv("HOME") ))
			fprintf( stderr, "%s: No such file or directory\n", s);
		}
		else {
		    while ( isspace( *s ) ) s++;
		    if (chdir( s ))
			fprintf( stderr, "%s: No such file or directory\n", s);
		}
	    }
	    /*free up the wasted space*/
	    if (par != NULL)
		free(par);
	    if (depth == -1)
		YYRETURN( back & PROMPT );
	}
#line 1639 "y.tab.c"
    break;

  case 14: /* line: endline  */
#line 383 "shell/parser.y"
        {
	    echo_bufptr = 0;

	    if (depth == -1)
		YYRETURN( (yyvsp[0].ival) );
	}
#line 1650 "y.tab.c"
    break;

  case 15: /* line: ENDFILE  */
#line 390 "shell/parser.y"
        {
	    echo_bufptr = 0;

	    if (depth == -1)
		YYRETURN(-1);
	}
#line 1661 "y.tab.c"
    break;

  case 16: /* line: QUIT  */
#line 397 "shell/parser.y"
        {
	    do_echo();

	    if (depth == -1)
		YYRETURN(-1);
	}
#line 1672 "y.tab.c"
    break;

  case 17: /* line: error endline  */
#line 404 "shell/parser.y"
        {
	    do_echo();

	    fprintf(stderr, "illegal input\n");
	    if (depth == -1)
		YYRETURN( (yyvsp[0].ival) & PROMPT );
	}
#line 1684 "y.tab.c"
    break;

  case 20: /* endline: EOL  */
#line 422 "shell/parser.y"
        {  (yyval.ival) = ~BACKGROUND & PROMPT; }
#line 1690 "y.tab.c"
    break;

  case 21: /* endline: BACK EOL  */
#line 424 "shell/parser.y"
        {  (yyval.ival) = BACKGROUND | PROMPT; }
#line 1696 "y.tab.c"
    break;

  case 22: /* endline: DELIMIT  */
#line 426 "shell/parser.y"
        {  (yyval.ival) = 0; }
#line 1702 "y.tab.c"
    break;

  case 23: /* endline: BACK DELIMIT  */
#line 428 "shell/parser.y"
        {  (yyval.ival) = BACKGROUND; }
#line 1708 "y.tab.c"
    break;

  case 24: /* rdir: REDIRECT  */
#line 433 "shell/parser.y"
        {
	    fprintf(stderr, "no file for redirection!\n");
	    (yyval.sval) = NULL;
	}
#line 1717 "y.tab.c"
    break;

  case 25: /* rdir: REDIRECT PARAMETER  */
#line 438 "shell/parser.y"
        {
	    (yyval.sval) = (yyvsp[0].sval);
	}
#line 1725 "y.tab.c"
    break;

  case 26: /* rdir: %empty  */
#line 442 "shell/parser.y"
        {
	    (yyval.sval) = NULL;
	}
#line 1733 "y.tab.c"
    break;

  case 27: /* parameter: PARAMETER  */
#line 449 "shell/parser.y"
        {
	    (yyval.sval) = (yyvsp[0].sval);
	}
#line 1741 "y.tab.c"
    break;

  case 28: /* parameter: %empty  */
#line 453 "shell/parser.y"
        {
	    (yyval.sval) = NULL;
	}
#line 1749 "y.tab.c"
    break;

  case 29: /* act1: %empty  */
#line 462 "shell/parser.y"
        {

	    /*increment the nesting depth for storage*/
	    depth++;

	    /*malloc off length if not done*/
	    if (store[depth] == NULL) {
		store_len[depth] = 1024;
		store[depth] = (char *)malloc(store_len[depth]);
		}
	    /*zero pointer*/
	    store_ptr[depth] = 0;

	    tp[depth] = prompt;
	    prompt = "> ";
	    te[depth] = echoall;
	    echoall = FALSE;

	}
#line 1773 "y.tab.c"
    break;

  case 30: /* act2: %empty  */
#line 486 "shell/parser.y"
        {
	    char *s;

	    *( store[depth] + store_ptr[depth]) = '\0';

	    s = (char *)malloc( strlen( store[depth] ) + 1 );
	    strcpy(s, store[depth]);

	    (yyval.sval) = s;
	}
#line 1788 "y.tab.c"
    break;


#line 1792 "y.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 499 "shell/parser.y"



int yyerror()
{ return 1; }
