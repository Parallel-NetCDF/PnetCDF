#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void main(){
	uint64_t a,b,c;
	b=1;
	for (a=1; a<64; a=2*2) 
		b = b*2;
	printf("a:%ld, b:%ld\n", a,b);

}
