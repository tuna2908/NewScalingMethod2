//
///* include -------------------------------------
//------------------------------------------------*/
//#include <conio.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//
///* define, constant ----------------------------------------------
//------------------------------------------------------------------*/
//#define M_PI 3.14159265359
//#define QUEUE_SIZE 1
//
//typedef float SAMPLE;
//
//typedef struct SIGNAL {
//	SAMPLE *raw_signal;
//	int frame_length;		//so sample trong 1 frame
//	int step_lengh;		//do dai bc nhay
//	int num_frame;			//so frame trong 1 audio signal
//	int signal_length;
//}SIGNAL;
//
//typedef struct COMPLEX {
//	float img;
//	float real;
//	float magnitude;
//}COMPLEX;
//
//typedef struct hyper_vector {
//	int row;
//	int col;
//	int dim;
//	SAMPLE *data;
//}hyper_vector;
//
///*
//functions
//*/
//
//SAMPLE * read_audio_signal_from_file(char * path, int* size)
//{
//	int i, j;
//	FILE *fp = fopen(path, "r");
//
//	fscanf(fp, "%d", size);
//
//	SAMPLE *audio_signal = (SAMPLE *)malloc(sizeof(SAMPLE) * (*size));
//	SAMPLE tmp;
//
//	if (fp == NULL) {
//		fprintf(stderr, "file no exist!!! \n");
//		exit(1);
//	}
//	else {
//		for (i = 0; i < (*size); ++i) {
//			fscanf(fp, "%f", &tmp);
//			audio_signal[i] = tmp;
//		}
//	}
//	fclose(fp);
//	return audio_signal;
//}
//
//float * butterworth_bandpass_v2(int order, float * signal, int size, float sample_rate, float high_freq_cutoff, float low_freq_cutoff, float *A,
//	float *d1, float *d2, float *d3, float *d4, float *w0, float *w1, float *w2, float *w3, float *w4, float *x)
//{
//	float a = cos(M_PI * (high_freq_cutoff + low_freq_cutoff) / sample_rate) / cos(M_PI * (high_freq_cutoff - low_freq_cutoff) / sample_rate);
//	float a_2 = a * a;
//	float b = tan(M_PI * (high_freq_cutoff - low_freq_cutoff) / sample_rate);
//	float b_2 = b * b;
//	float r;
//	for (int i = 0; i < order; ++i) {
//		r = sin(M_PI * (2.0 * i + 1.0) / (4.0 * order));
//		sample_rate = b_2 + 2.0 * b * r + 1.0;
//		A[i] = b_2 / sample_rate;
//		d1[i] = 4.0 * a * (1.0 + b * r) / sample_rate;
//		d2[i] = 2.0 * (b_2 - 2.0 * a_2 - 1.0) / sample_rate;
//		d3[i] = 4.0 * a * (1.0 - b * r) / sample_rate;
//		d4[i] = -(b_2 - 2.0 * b * r + 1.0) / sample_rate;
//	}
//	for (int i = 0; i < size; ++i) {
//		for (int j = 0; j < order; ++j) {
//			if (j == 0) {
//				w0[j] = d1[j] * w1[j] + d2[j] * w2[j] + d3[j] * w3[j] + d4[j] * w4[j] + signal[i];
//			}
//			else {
//				w0[j] = d1[j] * w1[j] + d2[j] * w2[j] + d3[j] * w3[j] + d4[j] * w4[j] + x[i];
//			}
//			x[i] = A[j] * (w0[j] - 2.0 * w2[j] + w4[j]);
//			w4[j] = w3[j];
//			w3[j] = w2[j];
//			w2[j] = w1[j];
//			w1[j] = w0[j];
//		}
//	}
//	return x;
//}
//
//float Window(float x) {
//	return (1.0*(1 + cos(2 * M_PI*x))) / 2;
//}
//
//main() {
//	char *path = "./data/0_0.txt";
//	int Nin, Nout;
//	float period_ratio = 1.5;
//
//	int inptr = 0, outprt = 0, period_len, old_zero_crossing = 0;
//	SAMPLE *out, oldx;
//	SAMPLE *in = read_audio_signal_from_file(path, &Nin);
//
//	Nout = Nin;
//	out = (SAMPLE *)calloc(Nin, sizeof(SAMPLE));
//
//	SAMPLE *temp = (SAMPLE*)calloc(1, sizeof(SAMPLE));
//	int order = 2;
//	float *A = (float *)malloc(sizeof(float) * order);
//	float *d1 = (float *)malloc(sizeof(float) * order);
//	float *d2 = (float *)malloc(sizeof(float) * order);
//	float *d3 = (float *)malloc(sizeof(float) * order);
//	float *d4 = (float *)malloc(sizeof(float) * order);
//	float *x = (float *)malloc(sizeof(float) * QUEUE_SIZE);
//	float *w0 = (float *)calloc(order, sizeof(float));
//	float *w1 = (float *)calloc(order, sizeof(float));
//	float *w2 = (float *)calloc(order, sizeof(float));
//	float *w3 = (float *)calloc(order, sizeof(float));
//	float *w4 = (float *)calloc(order, sizeof(float));
//
//	int dem = 0;
//	for (inptr = 0; inptr < Nin; inptr++) {
//		oldx = x[0];
//		temp[0] = in[inptr];
//
//		x = butterworth_bandpass_v2(2, temp, QUEUE_SIZE, 16000, 4000, 500, A, d1, d2, d3, d4, w0, w1, w2, w3, w4, x);
//
//
//		////zero crossing
//		//if (oldx>0 && x[0] <= 0 && dem == 0)
//		//{
//		//	old_zero_crossing = inptr;
//		//	dem++;
//		//}
//		//if (oldx > 0 && x[0] <= 0 && dem > 0)
//
//		if (oldx > 0 && x[0] <= 0)
//		{
//			period_len = inptr - old_zero_crossing;
//			old_zero_crossing = inptr;
//		anchor:
//			if ((1.0*outprt / Nout) < (1.0*inptr / Nin)) {
//				outprt += period_ratio * period_len;
//
//				for (int n = -period_len; n < period_len; n++) {
//					out[n + outprt] += in[n + inptr] * Window((1.0*n / period_len));
//				}
//				goto anchor;
//			}
//			else {
//				continue;
//			}
//		}
//		else {
//			continue;
//		}
//	}
//
//	FILE *fout = fopen("output.txt", "w");
//	for (int i = 0; i < Nout; i++) {
//		printf("%f\n", out[i]);
//		fprintf(fout, "%f\n", out[i]);
//	}
//	fclose(fout);
//	_getch();
//}






//float *realloc_same_add(float *data, int length, int new_length) {
//	for (int i = 0; i < length; ++i) {
//		printf("%d : %f\n", i, data[i]);
//	}
//	float *x = (float *)malloc(sizeof(float) * length);
//	for (int i = 0; i < length; ++i) {
//		x[i] = data[i];
//	}
//	free(data);
//	data = (float *)malloc(sizeof(float) * new_length);
//	for (int i = 0; i < length; ++i) {
//		data[i] = x[i];
//	}
//	free(x);
//	return data;
//}


//void ahi(char *tpath) {
//	int  label_cur = 0, i = 0;
//	char *label = (char *)malloc(sizeof(char) * 2);
//
//	int size = 0,max_index;
//	const char *default_ext = ".txt";
//
//	size_t len_path = strlen(tpath);
//	char *path = (char *)malloc(sizeof(char) * (len_path + 11));
//	char *path_conf = (char *)malloc(sizeof(char) * (len_path + 10 + 11));
//
//	len_path = strlen(path);
//
//	strcpy(path, tpath);
//	strcat(path, "realData\\");
//
//	strcpy(path_conf, path);
//	strcat(path_conf, "config.txt");
//
//	int dem = 0;
//
//	//load number of classes
//	FILE *fconf = fopen(path_conf, "r");
//	fscanf(fconf, "%d", &max_index);
//	fclose(fconf);
//
//	while (label_cur < (max_index + 1))
//	{
//		if (label_cur < 10) {
//			label = (char *)malloc(sizeof(char));
//		}
//		else if (label_cur > 9 && label_cur < 100) {
//			label = (char *)malloc(sizeof(char) * 2);
//		}
//
//		sprintf(label, "%d", label_cur);
//		size_t len_path_tmp = strlen(path) + strlen(label) + 1;
//		char *temp = (char*)malloc(sizeof(char) * len_path_tmp);
//		strcpy(temp, path);
//		strcat(temp, label);
//		strcat(temp, "_");
//		while (true) {
//			char *path_file;
//			char *index;
//			if (i < 10) { index = (char *)malloc(sizeof(char)); }
//			else if (i >= 10 && i < 100) {
//				index = (char *)malloc(sizeof(char) * 2);
//			}
//			else {
//				index = (char *)malloc(sizeof(char) * 3);
//			}
//			sprintf(index, "%d", i);
//			size_t len = len_path_tmp + strlen(index) + 5;
//			path_file = (char *)malloc(len * sizeof(char));
//			strcpy(path_file, temp);
//			strcat(path_file, index);
//			strcat(path_file, default_ext);
//			printf("path : %s \n", path_file);
//			if (check_path(path_file)) {
//				free(path_file);
//				break;
//			}
//			dem++;
//		}
//		++label_cur;
//		i = 0;
//	}
//	fprintf(fconf, "%d", dem);
//	fclose(fconf);
//}