#ifndef _BASE64_H
#define _BASE64_H

int b64_encode (char *dest,
		const unsigned char *src,
		int len);
int b64_decode (char *dest,
		const char *src);

void b64_decode_mio ( char *dest,  char *src );

void b64_decode_mio_max ( char *dest,  char *src , int max_bytes );

#endif /* BASE64_H */

