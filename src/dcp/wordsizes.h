#ifndef WORDSIZES_H_INC
#define WORDSIZES_H_INC

typedef char  S8BIT;		/* signed 8 bits */
typedef short S16BIT;		/* signed 16 bits */
typedef int   S32BIT;		/* signed 32 bits */
typedef unsigned char  U8BIT;	/* unsigned 8 bits */
typedef unsigned short U16BIT;	/* unsigned 16 bits */
typedef unsigned int   U32BIT;	/* unsigned 32 bits */

/***************************************/
/* stuff for compatibility with wave.h */
/***************************************/
#define WORD U16BIT
#define DWORD U32BIT

/**************************************************/
/* stuff for compatibility with aiff.h and aifc.h */
/**************************************************/
/*
 * Mac Type		Equivalence
 * ==============		===========
 * char			S8BIT
 * unsigned char		U8BIT
 * short			S16BIT
 * unsigned short		U16BIT
 * long			S32BIT
 * unsigned long		U32BIT
 * extended		80-bit IEEE-754 float; hope we can avoid this
 * 			fake this with 64-bit float followed by 2 NULL bytes
 * pstring			Pascal string; hope we can avoid this
 * 			first byte is length, followed by chars
 * 			pad with NULL if necessary to get even number of bytes
 * ID			char[4]
 * OSType			char[4]
 */
typedef struct {
    char bytes[10];
}F80BIT;		/* room to hold an IEEE 754 float with 80 bits */

/*
 * build a simple test program
 */
#ifdef BUILD_MAIN
main()
{
    printf ("sizeof(S8BIT)  = %d %s\n", sizeof(S8BIT),
	    (sizeof(S8BIT)==1?"OK":"ERROR!"));
    printf ("sizeof(S16BIT) = %d %s\n", sizeof(S16BIT),
	    (sizeof(S16BIT)==2?"OK":"ERROR!"));
    printf ("sizeof(S32BIT) = %d %s\n", sizeof(S32BIT),
	    (sizeof(S32BIT)==4?"OK":"ERROR!"));
    printf ("sizeof(U8BIT)  = %d %s\n", sizeof(U8BIT),
	    (sizeof(U8BIT)==1?"OK":"ERROR!"));
    printf ("sizeof(U16BIT) = %d %s\n", sizeof(U16BIT),
	    (sizeof(U16BIT)==2?"OK":"ERROR!"));
    printf ("sizeof(U32BIT) = %d %s\n", sizeof(U32BIT),
	    (sizeof(U32BIT)==4?"OK":"ERROR!"));
}
#endif /* BUILD_MAIN */
#endif /* WORDSIZES_H_INC */
