#include "stdio.h"
#include "JPGDecoder.h"

void main(char* args, char** kwargs) 
{
	Image* img = DecodeJPEG("Images/seven.jpg");
	printf("Iamge shape: (%d, %d, %d)\n", img->height, img->width, img->channels);

	for (int i = 0; i < img->height; i++) {
		for (int j = 0; j < img->width; j++) {
			int n_spaces = 0;
			if (img->data[i][j][0] < 10) {
				n_spaces = 3;
			}
			else if (img->data[i][j][0] < 100) {
				n_spaces = 2;
			}
			else if (img->data[i][j][0] < 300) {
				n_spaces = 1;
			}
			printf("%d", img->data[i][j][0]);
			for (int n = 0; n < n_spaces; n++) {
				printf(" ");
			}
		}
		printf("\n");
	}
}