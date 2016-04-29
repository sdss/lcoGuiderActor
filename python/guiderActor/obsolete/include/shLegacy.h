#ifndef SHLEGACY_H
#define SHLEGACY_H


#define S16 short int		//replacing dervish def
#define U8 unsigned char
#define   SH_SUCCESS		0x8001c009
#define SHERRORBASE		0x4000
#define SH_GENERIC_ERROR		SHERRORBASE	/* a generic error */

typedef struct mask{
   int nrow;               /* number of rows in mask */
   int ncol;               /* number of columns in mask */
   unsigned char **rows;   /* pointer to pointers to rows */

} MASK;

typedef struct region{
   int		  nrow;		/* number of rows in image */
   int		  ncol;		/* number of columns in image */
   S16 	      **rows_s16;
} REGION;

void shError(char *fmt, ...);
void shDebug(int level, char *fmt, ...);
void *shMalloc(size_t size);
void shFree(void *ptr);

#endif
