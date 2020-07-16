#ifndef SWAP_H_INC
#define SWAP_H_INC

#include "wordsizes.h"

#define swap_s16bit(a) { \
	union S16BITbytesUnion { S16BIT w; char c[2]; } s16bit_bytes; \
	S8BIT tmp_byte; \
	\
	s16bit_bytes.w = *(a); \
        tmp_byte = s16bit_bytes.c[0]; \
        s16bit_bytes.c[0] = s16bit_bytes.c[1]; \
        s16bit_bytes.c[1] = tmp_byte; \
        *(a) = s16bit_bytes.w; \
	}
#define swap_s32bit(a) { \
	union S32BITbytesUnion { S32BIT w; char c[4]; } s32bit_bytes; \
	S8BIT tmp_byte; \
	\
	s32bit_bytes.w = *(a); \
        tmp_byte = s32bit_bytes.c[0]; \
        s32bit_bytes.c[0] = s32bit_bytes.c[3]; \
        s32bit_bytes.c[3] = tmp_byte; \
        tmp_byte = s32bit_bytes.c[1]; \
        s32bit_bytes.c[1] = s32bit_bytes.c[2]; \
        s32bit_bytes.c[2] = tmp_byte; \
        *(a) = s32bit_bytes.w; \
	}
#define swap_u16bit(a) { \
	union U16BITbytesUnion { U16BIT w; char c[2]; } u16bit_bytes; \
	U8BIT tmp_ubyte; \
	\
	u16bit_bytes.w = *(a); \
        tmp_ubyte = u16bit_bytes.c[0]; \
        u16bit_bytes.c[0] = u16bit_bytes.c[1]; \
        u16bit_bytes.c[1] = tmp_ubyte; \
        *(a) = u16bit_bytes.w; \
	}
#define swap_u32bit(a) { \
	union U32BITbytesUnion { U32BIT w; char c[4]; } u32bit_bytes; \
	U8BIT tmp_ubyte; \
	\
	u32bit_bytes.w = *(a); \
        tmp_ubyte = u32bit_bytes.c[0]; \
        u32bit_bytes.c[0] = u32bit_bytes.c[3]; \
        u32bit_bytes.c[3] = tmp_ubyte; \
        tmp_ubyte = u32bit_bytes.c[1]; \
        u32bit_bytes.c[1] = u32bit_bytes.c[2]; \
        u32bit_bytes.c[2] = tmp_ubyte; \
        *(a) = u32bit_bytes.w; \
	}
#define swap_f64bit(a) { \
	union F64BITbytesUnion { F64BIT f; char c[8]; } f64bit_bytes; \
	U8BIT tmp_ubyte; \
	\
	f64bit_bytes.f = *(a); \
        tmp_ubyte = f64bit_bytes.c[0]; \
        f64bit_bytes.c[0] = f64bit_bytes.c[7]; \
        f64bit_bytes.c[7] = tmp_ubyte; \
        tmp_ubyte = f64bit_bytes.c[1]; \
        f64bit_bytes.c[1] = f64bit_bytes.c[6]; \
        f64bit_bytes.c[6] = tmp_ubyte; \
        tmp_ubyte = f64bit_bytes.c[2]; \
        f64bit_bytes.c[2] = f64bit_bytes.c[5]; \
        f64bit_bytes.c[5] = tmp_ubyte; \
        tmp_ubyte = f64bit_bytes.c[3]; \
        f64bit_bytes.c[3] = f64bit_bytes.c[4]; \
        f64bit_bytes.c[4] = tmp_ubyte; \
        *(a) = f64bit_bytes.f; \
	}
#define swap_f80bit(a) { \
	union F80BITbytesUnion { F80BIT f; char c[8]; } f80bit_bytes; \
	U8BIT tmp_ubyte; \
	\
	f80bit_bytes.f = *(a); \
        tmp_ubyte = f80bit_bytes.c[0]; \
        f80bit_bytes.c[0] = f80bit_bytes.c[9]; \
        f80bit_bytes.c[9] = tmp_ubyte; \
        tmp_ubyte = f80bit_bytes.c[1]; \
        f80bit_bytes.c[1] = f80bit_bytes.c[8]; \
        f80bit_bytes.c[8] = tmp_ubyte; \
        tmp_ubyte = f80bit_bytes.c[2]; \
        f80bit_bytes.c[2] = f80bit_bytes.c[7]; \
        f80bit_bytes.c[7] = tmp_ubyte; \
        tmp_ubyte = f80bit_bytes.c[3]; \
        f80bit_bytes.c[3] = f80bit_bytes.c[6]; \
        f80bit_bytes.c[6] = tmp_ubyte; \
        tmp_ubyte = f80bit_bytes.c[4]; \
        f80bit_bytes.c[4] = f80bit_bytes.c[5]; \
        f80bit_bytes.c[5] = tmp_ubyte; \
        *(a) = f80bit_bytes.f; \
	}

/* for wave.h compatibility */
#define swap_word(a)  swap_u16bit(a)
#define swap_dword(a) swap_u32bit(a)

#endif /* SWAP_H_INC */
