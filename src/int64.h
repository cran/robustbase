#ifndef __GNUC__
# define __extension__
#endif

/* This should work on Win64, as long is 4 bytes but long long is 8 bytes. */
#if defined __GNUC__ && __GNUC__ >= 2
__extension__ typedef long long int int64_t;
#else
typedef long long int int64_t;
/* FIXME: need double when long long is not available ? */
#endif
