#include "JPGDecoder.h"

static const size_t ZigZag[64] = {
	0,  1,  8, 16,  9,  2,  3, 10,
	17, 24, 32, 25, 18, 11,  4,  5,
	12, 19,	26, 33, 40, 48, 41, 34,
	27, 20, 13,  6,  7, 14, 21, 28,
	35, 42, 49, 56, 57, 50, 43, 36,
	29, 22, 15, 23, 30, 37, 44, 51,
	58, 59, 52, 45, 38, 31, 39, 46,
	53, 60, 61, 54, 47, 55, 62, 63
};

static BinaryData* ReadJPEG(char* filename) 
{
	int fd = fopen(filename, "rb");
	if (fd == -1) {
		printf("Can't open the file");
	}

	fseek(fd, 0, SEEK_END);
	long fsize = ftell(fd);
	fseek(fd, 0, SEEK_SET);

	char* data = malloc(fsize + 1);

	fread(data, fsize, 1, fd);
	fclose(fd);

	data[fsize] = '\0';

	BinaryData* bData = malloc(sizeof(BinaryData*));
	bData->data = data;
	bData->length = fsize + 1;

	return bData;
}

static void DefineQuantizationTables(char** chunk, int chunk_size, JPEG** jpeg) 
{
	int table_num = 0;
	register unsigned int iter = 0;
	while (iter < chunk_size) {
		table_num = (*chunk)[iter];
		(*jpeg)->quant_num = table_num + 1;
		(*jpeg)->quant = (unsigned int**)realloc((*jpeg)->quant, (table_num + 1) * sizeof(unsigned int*));
		(*jpeg)->quant[table_num] = (unsigned int*)malloc(64 * sizeof(unsigned int));
		for (register size_t i = iter; i < iter + 64; i++) {
			(*jpeg)->quant[table_num][i - iter] = (*chunk)[i + 1];
		}
		iter += 65;
	}
	free((*chunk));
}

static void BaselineDCT(char** chunk, int chunk_size, JPEG** jpeg) 
{
	unsigned int hdr = (*chunk)[0];
	(*jpeg)->height = (unsigned char)(*chunk)[1] << 8 | (unsigned char)(*chunk)[2];
	(*jpeg)->width = (unsigned char)(*chunk)[3] << 8 | (unsigned char)(*chunk)[4];
	(*jpeg)->channels = (*chunk)[5];

	(*jpeg)->padded_height = (*jpeg)->height;
	(*jpeg)->padded_width = (*jpeg)->width;

	if ((*jpeg)->height % 8 != 0) {
		(*jpeg)->padded_height = ((unsigned int)((*jpeg)->height / 8) + 1) * 8;
	}
	if ((*jpeg)->width % 8 != 0) {
		(*jpeg)->padded_width = ((unsigned int)((*jpeg)->width / 8) + 1) * 8;
	}

	(*jpeg)->image = (int*)calloc((*jpeg)->padded_height * (*jpeg)->padded_width, sizeof(int));
	(*jpeg)->h_sparse_coefs = (unsigned int*)calloc((*jpeg)->channels, sizeof(unsigned int));
	(*jpeg)->v_sparse_coefs = (unsigned int*)calloc((*jpeg)->channels, sizeof(unsigned int));
	(*jpeg)->quantMapping = (unsigned int*)calloc((*jpeg)->channels, sizeof(unsigned int));

	(*jpeg)->H_MAX = 0;
	(*jpeg)->V_MAX = 0;

	for (register size_t i = 0; i < (*jpeg)->channels; i++) {
		unsigned int id = (*chunk)[6 + i * 3];
		unsigned int sparse = (*chunk)[7 + i * 3];
		(*jpeg)->h_sparse_coefs[i] = (sparse & 0xf0) >> 4;
		(*jpeg)->v_sparse_coefs[i] = sparse & 0x0f;
		(*jpeg)->H_MAX = max((*jpeg)->H_MAX, (*jpeg)->h_sparse_coefs[i]);
		(*jpeg)->V_MAX = max((*jpeg)->V_MAX, (*jpeg)->h_sparse_coefs[i]);
		(*jpeg)->quantMapping[i] = (*chunk)[8 + i * 3];
	}

	free((*chunk));
}

static void DefineHuffmanTables(char** chunk, int chunk_size, JPEG** jpeg) 
{
	unsigned int off = 0;
	while (off < chunk_size) {
		unsigned int hdr = (*chunk)[off];
		unsigned int DC_AC = (hdr & 0xf0) >> 4;
		unsigned int tId = hdr & 0x0f;
		off += 1;

		unsigned int* lengths = (unsigned int*)malloc(16 * sizeof(unsigned int));
		unsigned int elem_num = 0;
		for (register size_t i = 0; i < 16; i++) {
			lengths[i] = (*chunk)[off + i];
			elem_num += lengths[i];
		}

		off += 16;
		unsigned int* elements = (unsigned int*)malloc(elem_num * sizeof(unsigned int));
		register size_t ii = 0;
		register size_t i = 0;
		while (ii < elem_num) {
			for (register size_t j = 0; j < lengths[i]; j++) {
				elements[ii + j] = (unsigned char)(*chunk)[off + j];
			}
			ii += lengths[i];
			off += lengths[i];
			i += 1;
		}

		if (!DC_AC) {
			(*jpeg)->ht_num += 1;
			(*jpeg)->huffman_tables_dc = (HuffManTable**)realloc((*jpeg)->huffman_tables_dc, (tId + 1) * sizeof(HuffManTable*));
		}
		else {
			(*jpeg)->huffman_tables_ac = (HuffManTable**)realloc((*jpeg)->huffman_tables_ac, (tId + 1) * sizeof(HuffManTable*));
		}
		HuffManTable* hf = malloc(sizeof(HuffManTable));
		hf->DC_AC = DC_AC;
		hf->tId = tId;
		HuffManRoot* head = malloc(sizeof(HuffManRoot));
		head->length = 0;
		head->left_end = 0;
		head->right_end = 0;
		hf->head = head;
		GetHuffManBits(hf, lengths, elements);
		if (!DC_AC) {
			(*jpeg)->huffman_tables_dc[tId] = hf;
		}
		else {
			(*jpeg)->huffman_tables_ac[tId] = hf;
		}
	}

	free(*chunk);
}

static void GetHuffManBits(HuffManTable* hf, unsigned int* lengths, unsigned int* elements) 
{
	register size_t ii = 0;
	for (register size_t i = 0; i < 16; i++) {
		for (register size_t j = 0; j < lengths[i]; j++) {
			BitsFromLengths(hf->head, elements[ii], i, 0);
			ii += 1;
		}
	}
}

static unsigned int BitsFromLengths(HuffManRoot* root, unsigned int element, unsigned int pos, unsigned int end) 
{
	if (end == 0) {
		if (pos == 0) {
			if (root->length < 1) {
				root->left = element;
				root->left_end = 1;
				root->length += 1;
				return 1;
			}
			else if (root->length < 2) {
				root->right = element;
				root->right_end = 1;
				root->length += 1;
				return 1;
			}

			return 0;
		}

		if (root->length == 0) {
			HuffManRoot* new_root = malloc(sizeof(HuffManRoot));
			new_root->length = 0;
			new_root->left_end = 0;
			new_root->right_end = 0;
			root->left = new_root;
			root->length += 1;
		}

		if (BitsFromLengths(root->left, element, pos - 1, root->left_end) == 1) {
			return 1;
		}

		if (root->length == 1) {
			HuffManRoot* new_root = malloc(sizeof(HuffManRoot));
			new_root->length = 0;
			new_root->left_end = 0;
			new_root->right_end = 0;
			root->right = new_root;
			root->length += 1;
		}

		if (BitsFromLengths(root->right, element, pos - 1, root->right_end) == 1) {
			return 1;
		}
	}
	return 0;
}

static void freeHuffmanRoot(HuffManRoot* node) 
{
	if ((node->left_end) | (node->right_end)) {
		return;
	}

	else {
		freeHuffmanRoot(node->left);
		freeHuffmanRoot(node->right);
		free(node);
	}
}

static void freeHuffmanTable(HuffManTable* table) 
{
	freeHuffmanRoot(table->head);
	free(table);
}

static void BuildMatrix(double* mat, Stream* st, JPEG** jpeg, unsigned int idx, unsigned int* quant, int* dccoef) 
{
	unsigned int code = GetCode(st, (*jpeg)->huffman_tables_dc[idx]);
	unsigned int bits = GetBitN(st, code);
	*dccoef += DecodeNumber(code, bits);

	AddZigZag(mat, 0, (*dccoef) * quant[0]);

	register size_t l = 1;
	while (l < 64) {
		unsigned int code = GetCode(st, (*jpeg)->huffman_tables_ac[idx]);
		if (code == 0) {
			break;
		}
		if (code > 15) {
			l += (code >> 4);
			code &= 0xf;
		}

		unsigned int bits = GetBitN(st, code);

		if (l < 64) {
			int coef = DecodeNumber(code, bits);
			AddZigZag(mat, l, coef * quant[l]);
			l += 1;
		}
	}

	return mat;
}

static void BuildMCU(double* MCU, Stream* st, JPEG** jpeg, unsigned int idx, unsigned int* quant, int* dccoef, int cn) 
{
	for (register int i = 0; i < (*jpeg)->h_sparse_coefs[cn]; i++) {
		double mat[64] = { 0.0 };
		BuildMatrix(mat, st, jpeg, idx, quant, dccoef);

		for (register int j = i * 8 * (*jpeg)->h_sparse_coefs[cn];
			j < 8 * (*jpeg)->v_sparse_coefs[cn] + i * 8 * (*jpeg)->h_sparse_coefs[cn];
			j += (*jpeg)->v_sparse_coefs[cn]) 
		{
			for (register int k = 0; k < 8; k++) {
				MCU[j * 8 + k] = mat[(int)(j/(*jpeg)->v_sparse_coefs[cn]) * 8 - i * 64 + k];
			}
		}
	}

	if ((*jpeg)->v_sparse_coefs[cn] > 1) {
		for (register int i = 0; i < (*jpeg)->v_sparse_coefs[cn]; i++) {
			double mat[64] = { 0.0 };
			BuildMatrix(mat, st, jpeg, idx, quant, dccoef);

			for (register int j = i * 8 * (*jpeg)->h_sparse_coefs[cn] + 1;
				j < 8 * (*jpeg)->v_sparse_coefs[cn] + i * 8 * (*jpeg)->h_sparse_coefs[cn] + 1;
				j += (*jpeg)->v_sparse_coefs[cn])
			{
				for (register int k = 0; k < 8; k++) {
					MCU[j * 8 + k] = mat[(int)(j / (*jpeg)->v_sparse_coefs[cn]) * 8 - i * 64 + k];
				}
			}
		}
	}
}

static unsigned int StartOfScan(BinaryData* data, int hdrlen, JPEG** jpeg) 
{
	int new_size = data->length - hdrlen;
	char* slice = malloc(new_size);
	memcpy(slice, data->data + hdrlen, new_size);
	free(data->data);
	data->data = slice;
	data->length = new_size;
	RemoveFF00(data);

	Stream* stream = (Stream*)malloc(sizeof(Stream));
	stream->data = data->data;
	stream->pos = 0;

	(*jpeg)->image = (unsigned char*)calloc((*jpeg)->padded_height * (*jpeg)->padded_width * (*jpeg)->channels, sizeof(unsigned char));

	if ((*jpeg)->channels == 1) {
		int dc_coef = 0;
		for (int y = 0; y < (*jpeg)->padded_height / (8 * (*jpeg)->V_MAX); y++) {
			for (int x = 0; x < (*jpeg)->padded_width / (8 * (*jpeg)->H_MAX); x++) {
				double* MCU = (double*)calloc(64 * (*jpeg)->h_sparse_coefs[0] * (*jpeg)->v_sparse_coefs[0], sizeof(double));

				BuildMCU(MCU, stream, jpeg, 0, (*jpeg)->quant[(*jpeg)->quantMapping[0]], &dc_coef, 0);

				for (register size_t yy = 0; yy < 8 * (*jpeg)->V_MAX; yy++) {
					for (register size_t xx = 0; xx < 8 * (*jpeg)->H_MAX; xx++) {
						(*jpeg)->image[(x * 8 * (*jpeg)->H_MAX + xx) + ((y * 8 * (*jpeg)->V_MAX + yy) * (*jpeg)->padded_width)] = 
							GrayScaleConversion(MCU[(int)(yy / ((*jpeg)->H_MAX / (*jpeg)->h_sparse_coefs[0])) + 
								(int)(xx / ((*jpeg)->V_MAX / (*jpeg)->v_sparse_coefs[0]))*8*(*jpeg)->h_sparse_coefs[0]]);
					}
				}
				free(MCU);
			}
		}
	}

	else if ((*jpeg)->channels == 3) {
		int lumbdccoef = 0;
		int Cbdccoef = 0;
		int Crdccoef = 0;
		for (int y = 0; y < (*jpeg)->padded_height / (8 * (*jpeg)->V_MAX); y++) {
			for (int x = 0; x < (*jpeg)->padded_width / (8 * (*jpeg)->H_MAX); x++) {
				double* matL = (double*)calloc(64 * (*jpeg)->h_sparse_coefs[0] * (*jpeg)->v_sparse_coefs[0], sizeof(double));
				double* matCr = (double*)calloc(64 * (*jpeg)->h_sparse_coefs[1] * (*jpeg)->v_sparse_coefs[1], sizeof(double));
				double* matCb = (double*)calloc(64 * (*jpeg)->h_sparse_coefs[2] * (*jpeg)->v_sparse_coefs[2], sizeof(double));

				BuildMCU(matL, stream, jpeg, 0, (*jpeg)->quant[(*jpeg)->quantMapping[0]], &lumbdccoef, 0);
				BuildMCU(matCr, stream, jpeg, 1, (*jpeg)->quant[(*jpeg)->quantMapping[1]], &Crdccoef, 1);
				BuildMCU(matCb, stream, jpeg, 1, (*jpeg)->quant[(*jpeg)->quantMapping[2]], &Cbdccoef, 2);

				for (register size_t yy = 0; yy < 8 * (*jpeg)->V_MAX; yy++) {
					for (register size_t xx = 0; xx < 8 * (*jpeg)->H_MAX; xx++) {
						unsigned char pixel[3];
						ColorConversion(matL[(int)(yy / ((*jpeg)->H_MAX / (*jpeg)->h_sparse_coefs[0])) +
							(int)(xx / ((*jpeg)->V_MAX / (*jpeg)->v_sparse_coefs[0])) * 8 * (*jpeg)->h_sparse_coefs[0]],
							matCb[(int)(yy / ((*jpeg)->H_MAX / (*jpeg)->h_sparse_coefs[1])) +
							(int)(xx / ((*jpeg)->V_MAX / (*jpeg)->v_sparse_coefs[1])) * 8 * (*jpeg)->h_sparse_coefs[1]],
							matCr[(int)(yy / ((*jpeg)->H_MAX / (*jpeg)->h_sparse_coefs[2])) +
							(int)(xx / ((*jpeg)->V_MAX / (*jpeg)->v_sparse_coefs[2])) * 8 * (*jpeg)->h_sparse_coefs[2]],
							pixel);

						for (register size_t zz = 0; zz < 3; zz++) {
							(*jpeg)->image[(x * 8 * (*jpeg)->H_MAX + xx) + ((y * 8 * (*jpeg)->V_MAX + yy) * (*jpeg)->padded_width) 
								+ zz] = pixel[zz];
						}
					}
				}

				free(matL);
				free(matCb);
				free(matCr);
			}
		}
	}

	free(stream);

	return data->length;
}

static void RemoveFF00(BinaryData* data) 
{
	char* datapro = (char*)malloc(data->length);
	register unsigned int i = 0;
	register unsigned int j = 0;

	while (i < data->length) {
		unsigned char b = (unsigned char)data->data[i];
		unsigned char b_next = (unsigned char)data->data[i + 1];

		if (b == 0xff) {
			if (b_next != 0) {
				break;
			}
			datapro[j] = data->data[i];
			j += 1;
			i += 2;
		}

		else {
			datapro[j] = data->data[i];
			j += 1;
			i += 1;
		}
	}

	int new_size = j;
	free(data->data);
	data->data = (char*)malloc(new_size);
	memcpy(data->data, datapro, new_size);
	free(datapro);
	data->length = new_size;
}

static unsigned int GetCode(Stream* st, HuffManTable* ht) 
{
	unsigned int res = -1;
	while (1) {
		HuffManRoot* root = ht->head;
		while (1) {
			unsigned int bit = GetBit(st);
			if (bit) {
				if (root->right_end) {
					res = (unsigned int)root->right;
					break;
				}
				root = (HuffManRoot*)root->right;
			}

			else {
				if (root->left_end) {
					res = (unsigned int)root->left;
					break;
				}
				root = (HuffManRoot*)root->left;
			}
		}
		if (res != -1) {
			return res;
		}
	}
}

static unsigned int GetBit(Stream* st) 
{
	unsigned int b = st->data[st->pos >> 3];
	unsigned int s = 7 - (st->pos & 0x7);
	st->pos += 1;

	return (b >> s) & 1;
}

static unsigned int GetBitN(Stream* st, unsigned int code) 
{
	unsigned int val = 0;
	for (int i = 0; i < code; i++) {
		val = val * 2 + GetBit(st);
	}

	return val;
}

static int DecodeNumber(unsigned int code, unsigned int bits) 
{
	int l = pow(2, code-1);
	if (bits >= l) {
		return bits;
	}
	else {
		int d = pow(2, code) - 1;
		return bits - d;
	}
}

static void AddZigZag(double* mat, unsigned int zi, int coef) 
{
	unsigned int i = ZigZag[zi];
	unsigned int n = i & 0x7;
	unsigned int m = i >> 3;

	double an = NormCoef(n);
	double am = NormCoef(m);

	for (register size_t y = 0; y < 8; y++) {
		for (register size_t x = 0; x < 8; x++) {
			double nn = an * cos(n * M_PI * (x + 0.5) / 8.0);
			double mm = am * cos(m * M_PI * (y + 0.5) / 8.0);
			mat[y + x*8] += (nn*mm*coef);
		}
	}
}

static double NormCoef(unsigned int n) 
{
	return n == 0 ? sqrt(1.0 / 8.0) : sqrt(2.0 / 8.0);
}

static unsigned char Clamp(double val) 
{
	if (val > 255) {
		val = 255;
	}

	if (val < 0) {
		val = 0;
	}

	return (unsigned char)val;
}

static unsigned char GrayScaleConversion(double val) 
{
	return Clamp(val + 128);
}

static void ColorConversion(double Y, double Cr, double Cb, unsigned char* pixel) 
{
	double R = Cr * (2 - 2 * 0.299) + Y;
	double B = Cb * (2 - 2 * 0.114) + Y;
	double G = (Y - 0.114 * B - 0.299 * R) / 0.587;

	pixel[0] = Clamp(R + 128);
	pixel[1] = Clamp(B + 128);
	pixel[2] = Clamp(G + 128);
}

Image* DecodeJPEG(char* filename) 
{
	BinaryData* bData = ReadJPEG(filename);
	JPEG* jpeg = (JPEG*)malloc(sizeof(JPEG));
	jpeg->quant = (unsigned int**)malloc(sizeof(unsigned int*));
	jpeg->huffman_tables_dc = (HuffManTable**)malloc(sizeof(HuffManTable*));
	jpeg->huffman_tables_ac = (HuffManTable**)malloc(sizeof(HuffManTable*));
	jpeg->ht_num = 0;

	while (1) {
		unsigned char b1 = bData->data[0];
		unsigned char b2 = bData->data[1];

		unsigned int hdr = (b1 << 8 | b2);

		int lenchunk;

		if (hdr == 0xffd8) {
			lenchunk = 2;
		}

		else if (hdr == 0xffd9) {
			break;
		}

		else {
			unsigned char b1 = bData->data[2];
			unsigned char b2 = bData->data[3];

			lenchunk = (b1 << 8 | b2);
			lenchunk += 2;

			int chunk_size = lenchunk - 4;
			char* chunk = malloc(chunk_size);
			memcpy(chunk, bData->data + 4, chunk_size);
			if (hdr == 0xffdb) {
				DefineQuantizationTables(&chunk, chunk_size, &jpeg);
			}
			else if (hdr == 0xffc0) {
				BaselineDCT(&chunk, chunk_size, &jpeg);
			}
			else if (hdr == 0xffc4) {
				DefineHuffmanTables(&chunk, chunk_size, &jpeg);
			}
			else if (hdr == 0xffda) {
				lenchunk = StartOfScan(bData, lenchunk, &jpeg);
			}
		}

		int new_size = bData->length - lenchunk;
		if (new_size == 0) {
			free(bData->data);
			Image* img = (Image*)malloc(sizeof(Image));
			img->height = jpeg->height;
			img->width = jpeg->width;
			img->channels = jpeg->channels;

			img->data = (unsigned char***)malloc(img->height * sizeof(unsigned char**));
			for (register int i = 0; i < img->height; i++) {
				img->data[i] = (unsigned char**)malloc(img->width * sizeof(unsigned char*));
				for (register int j = 0; j < img->width; j++) {
					img->data[i][j] = (unsigned char*)malloc(img->channels * sizeof(unsigned char));
					for (register int k = 0; k < img->channels; k++) {
						img->data[i][j][k] = jpeg->image[i * jpeg->padded_width + j + k];
					}
				}
			}

			free(jpeg->image);
			for (register size_t i = 0; i < jpeg->quant_num; i++) {
				free(jpeg->quant[i]);
			}
			free(jpeg->quant);
			free(jpeg->quantMapping);

			free(jpeg->h_sparse_coefs);
			free(jpeg->v_sparse_coefs);

			for (register size_t i = 0; i < jpeg->ht_num; i++) {
				freeHuffmanTable(jpeg->huffman_tables_dc[i]);
				freeHuffmanTable(jpeg->huffman_tables_ac[i]);
			}
			free(jpeg->huffman_tables_dc);
			free(jpeg->huffman_tables_ac);
			free(jpeg);

			return img;
		}
		char* slice = malloc(new_size);
		memcpy(slice, bData->data + lenchunk, new_size);
		free(bData->data);
		bData->data = slice;
		bData->length = new_size;
	}
}