#include "stdio.h"
#include "JPGDecoder.h"

void main(char* args, char** kwargs) 
{
	Image* img = DecodeJPEG("Images/cat.jpg");
	printf("Iamge shape: (%d, %d, %d)", img->height, img->width, img->channels);
}