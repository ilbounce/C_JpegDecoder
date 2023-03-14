#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define M_PI 3.14159265358979323846

typedef struct {
	char* data;
	long length;
} BinaryData;

typedef struct {
	void* left;
	void* right;
	unsigned int left_end;
	unsigned int right_end;
	unsigned int length;
} HuffManRoot;

typedef struct {
	unsigned int DC_AC;
	unsigned int tId;
	HuffManRoot* head;
} HuffManTable;

typedef struct {
	char* data;
	int pos;
} Stream;

typedef struct {
	unsigned int height;
	unsigned int width;
	unsigned int channels;
	unsigned int padded_height;
	unsigned int padded_width;
	unsigned int** quant;
	unsigned int quant_num;
	unsigned int H_MAX;
	unsigned int V_MAX;
	unsigned int* h_sparse_coefs;
	unsigned int* v_sparse_coefs;
	unsigned int* quantMapping;
	unsigned char* image;
	HuffManTable** huffman_tables_dc;
	HuffManTable** huffman_tables_ac;
	unsigned int ht_num;
} JPEG;

typedef struct {

	unsigned int height;
	unsigned int width;
	unsigned int channels;
	unsigned char*** data;

} Image;

Image* DecodeJPEG(char* filename);

static BinaryData* ReadJPEG(char* filename);

static void DefineQuantizationTables(char** chunk, int chunk_size, JPEG** jpeg);

static void BaselineDCT(char** chunk, int chunk_size, JPEG** jpeg);

static void DefineHuffmanTables(char** chunk, int chunk_size, JPEG** jpeg);

static void GetHuffManBits(HuffManTable* hf, unsigned int* lengths, unsigned int* elements);

static unsigned int BitsFromLengths(HuffManRoot* root, unsigned int element, unsigned int pos, unsigned int end);

static void freeHuffmanRoot(HuffManRoot* node);

static void freeHuffmanTable(HuffManTable* table);

static unsigned int StartOfScan(BinaryData* data, int hdrlen, JPEG** jpeg);

static void RemoveFF00(BinaryData* data);

static unsigned int GetBit(Stream* st);

static unsigned int GetBitN(Stream* st, unsigned int code);

static unsigned int GetCode(Stream* st, HuffManTable* ht);

static int DecodeNumber(unsigned int code, unsigned int bits);

static void AddZigZag(double* mat, unsigned int zi, int coef);

static double NormCoef(unsigned int n);

static void BuildMatrix(double* mat, Stream* st, JPEG** jpeg, unsigned int idx, unsigned int* quant, int* dccoef);

static void BuildMCU(double* MCU, Stream* st, JPEG** jpeg, unsigned int idx, unsigned int* quant, int* dccoef, int cn);

static unsigned char GrayScaleConversion(double val);

static unsigned char Clamp(double val);

static void ColorConversion(double Y, double Cr, double Cb, unsigned char* pixel);