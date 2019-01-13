/**************  Effects Program  *********************/

#include "LentPitShift.h"
#include "kiss_fft.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include<conio.h>
using namespace stk;

#define PI 3.14159265358979323846
#define TWOPI PI*2


void mainProcess(StkFloat shift_period, int numData, char* path);
void time_stretch(StkFloat *signal, int winsize, int fftsize, int shop, int ahop, float *CENTERFREQ, kiss_fft_cpx *Y, kiss_fft_cpx * cx_out_i,
	kiss_fft_cfg cfg_i, kiss_fft_cpx *cx_in, kiss_fft_cpx *cx_out, kiss_fft_cfg cfg, float *framed, float *Mag, float *Pha,
	float *PhaSy, float *old_Pha, float *old_PhaSy, float *PhaseDiff, float *dphi, float *freq, float *Y_out,
	int size, StkFloat *Out, int OutLen, int num_win);

float HanningWindow(float a, int frameLength)
{

	return 0.5*(1 - cos(2 * PI * a / (frameLength - 1)));
}

float HammingWindow(float a, int frameLength)
{
	return 0.54 - 0.46 * cos(2 * PI * a / (frameLength - 1));
}

float magnitude(float real, float img)
{
	return sqrt(real*real + img * img);
}


int check_path(char * path)
{
	FILE *fp = fopen(path, "r");
	if (fp == NULL) {
		fprintf(stderr, "file no exist!!! \n");
		return 1;
	}
	fclose(fp);
	return 0;
}

void usage() {
	printf(" executable.exe -shift (recommend <0.2) number_of_scaleData (ex: 100) .\\folderName\\");
}

int find_args(int argc, char * argv[], char * arg)
{
	int i;
	for (i = 0; i < argc; ++i) {
		if (!argv[i]) continue;
		if (0 == strcmp(argv[i], arg)) {
			return 1;
		}
	}
	return 0;
}

int valid(char *path, char * num) {
	if (path == NULL || num == NULL) {
		printf("DD");
		return 0;
	}
	return 1;
}

void test_r() {
	StkFloat tempf = 0, step, shift_period = 0.2;
	StkFloat from = 1 + shift_period, to = 1 - shift_period;
	time_t t;
	//int size;
	step = (from - to) / 500;

	srand((unsigned)time(&t));
	for (int n = 0; n < 500; n++) {
		printf("%lf\n", (from - (rand() % (500 + 1))*step));
	}

	_getch();
}

void test_path() {
	int len = strlen("D:\\Work_place_game\\data\\");
	char *path = (char*)malloc(sizeof(char)*len);
	strcpy(path, "D:\\Work_place_game\\data\\");
	mainProcess(0.2, 100, path);
	_getch();
}

void test_pitch() {
	/*LentPitShift lshifter;
	StkFloat tempf = 0, step;
	int size;
	lshifter.setPeriod(0.7);
	fscanf(fin, "%d", &size);
	for (int i = 0; i < size; i++) {
		fscanf(fin, "%lf", &tempf);
		lshifter.tick(tempf, fout);
	}
	fclose(fin);
	fclose(fout);
	_getch();*/

	LentPitShift lshifter;
	StkFloat temp;
	int size;
	FILE *fin = fopen("D:\\Work_place_game\\data\\realData\\0_0.txt", "r");
	FILE *fout = fopen("D:\\Work_place_game\\data\\scaleData\\0_0.txt", "w");

	fscanf(fin, "%d", &size);
	lshifter.setPeriod(0.75);
	for (int i = 0; i < size; i++) {
		fscanf(fin, "%lf", &temp);
		StkFloat sample = lshifter.tick(temp, fout);
		//fprintf(fout,"%.9lf\n", sample);

		//printf("%lf - %.9lf\n", temp, sample);


	}
	fclose(fout);
	fclose(fin);
}


int main(int argc, char *argv[])
{
	//Variables
	//test_pitch();
	//_getch();
	//return 1;
	//test_path();
	//return 1;
	char *path = argv[argc - 1];
	char *numDataString = argv[argc - 2];
	int numData = atoi(numDataString);

	//int shifted = find_args(argc, argv, "-shift");
	StkFloat shift_period = 0.2;

	if (valid(path, numDataString))
	{
		for (int i = 1; i < argc; i++) {
			if (!strcmp(argv[i], "-shift")) {
				shift_period = atoi(argv[++i]);
				break;
			}
		}

		//printf("%d - %s", numData, path);
		//Process
		//_getch();
		mainProcess(shift_period, numData, path);
		return 1;
	}
	usage();
}

void mainProcess(StkFloat shift_period, int numData, char* tpath) {
	LentPitShift lshifter;
	StkFloat  step;
	StkFloat from = 1 + shift_period, to = 1 - shift_period;

	time_t t;
	//int size;
	step = (from - to) / 300;

	/////////TIME STRETCHTER PARAMS////////////////////////////
	int winsize = 512;
	int fftsize = winsize;
	int shop = winsize / 4;
	float *CENTERFREQ = (float*)malloc(sizeof(float)*fftsize);
	kiss_fft_cpx *Y = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*fftsize);
	kiss_fft_cpx * cx_out_i = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*fftsize);
	kiss_fft_cfg cfg_i = kiss_fft_alloc(fftsize, 1, 0, 0);
	kiss_fft_cpx * cx_in = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*fftsize);
	kiss_fft_cpx * cx_out = (kiss_fft_cpx*)malloc(sizeof(kiss_fft_cpx)*fftsize);
	kiss_fft_cfg cfg = kiss_fft_alloc(fftsize, 0, 0, 0);
	float *framed = (float*)malloc(sizeof(float)*winsize);
	float *Mag = (float*)malloc(sizeof(float)*fftsize);
	float *Pha = (float*)malloc(sizeof(float)*fftsize);
	float *PhaSy = (float*)malloc(sizeof(float)*fftsize);
	float *old_Pha = (float*)malloc(sizeof(float)*fftsize);
	float *old_PhaSy = (float*)malloc(sizeof(float)*fftsize);
	float *PhaseDiff = (float*)malloc(sizeof(float)*fftsize);
	float *dphi = (float*)malloc(sizeof(float)*fftsize);
	float *freq = (float*)malloc(sizeof(float)*fftsize);
	float *Y_out = (float*)malloc(sizeof(float)*fftsize);

	///////////////////////////////////////////////////////////////////

	int  label_cur = 0, i = 0, bound = 0;
	char *label = (char *)malloc(sizeof(char) * 2);

	int size = 0, max_index, sum = 0;
	const char *default_ext = ".txt";

	size_t len_path = strlen(tpath);
	char *path = (char *)malloc(sizeof(char) * (len_path + 11));
	char *path2 = (char *)malloc(sizeof(char) * (len_path + 11));

	char *path_conf = (char *)malloc(sizeof(char) * (len_path + 10 + 11));
	char *path_info = (char *)malloc(sizeof(char) * (len_path + 10 + 11));



	len_path = strlen(path);

	strcpy(path, tpath);
	strcat(path, "realData\\");

	strcpy(path2, tpath);
	strcat(path2, "scaleData\\");

	strcpy(path_conf, path);
	strcat(path_conf, "config.txt");

	strcpy(path_info, path2);
	strcat(path_info, "info.txt");

	int dem = 0;

	//load number of classes
	FILE *fconf = fopen(path_conf, "r");
	fscanf(fconf, "%d", &max_index);
	fclose(fconf);

	FILE *finf = fopen(path_info, "w");


	while (label_cur < (max_index + 1))
	{
		if (label_cur < 10) {
			label = (char *)malloc(sizeof(char));
		}
		else if (label_cur > 9 && label_cur < 100) {
			label = (char *)malloc(sizeof(char) * 2);
		}

		sprintf(label, "%d", label_cur);
		size_t len_path_tmp = strlen(path) + strlen(label) + 1;
		char *temp = (char*)malloc(sizeof(char) * len_path_tmp);
		strcpy(temp, path);
		strcat(temp, label);
		strcat(temp, "_");
		while (true) {
			char *path_file;
			char *index;
			if (i < 10) { index = (char *)malloc(sizeof(char)); }
			else if (i >= 10 && i < 100) {
				index = (char *)malloc(sizeof(char) * 2);
			}
			else {
				index = (char *)malloc(sizeof(char) * 3);
			}
			sprintf(index, "%d", i);
			size_t len = len_path_tmp + strlen(index) + 5;
			path_file = (char *)malloc(len * sizeof(char));
			strcpy(path_file, temp);
			strcat(path_file, index);
			strcat(path_file, default_ext);
			printf("path : %s \n", path_file);
			if (check_path(path_file)) {
				free(path_file);
				sum += bound;
				fprintf(finf, "%d ", bound);
				bound = 0;
				break;
			}
			else {
				FILE * fout;
				FILE * fin = fopen(path_file, "r");
				fscanf(fin, "%d", &size);
				StkFloat *tempf = (StkFloat*)malloc(sizeof(StkFloat)*size);
				for (int i = 0; i < size; i++) {
					fscanf(fin, "%lf", &tempf[i]);
				}
				fclose(fin);
				int random = 0;
				int new_size;
				float factor;
				int ahop;
				int num_win;
				int SignalLen = size;

				char *path_file2;
				char *index2 = (char*)malloc(sizeof(char) * 3);
				srand((unsigned)time(&t));
				for (int n = bound; n < numData + bound; n++) {
					////TIME STRETCHING////
					random = rand() % (300 + 1);
					factor = from - random * step;
					ahop = floor(shop / factor);
					num_win = floor((SignalLen - winsize) / ahop);

					int OutLen = (num_win - 1)*shop + winsize;
					StkFloat *Out= (StkFloat*)calloc(OutLen, sizeof(StkFloat));

					time_stretch(tempf, winsize, fftsize, shop, ahop, CENTERFREQ, Y, cx_out_i, cfg_i, cx_in, cx_out, cfg, framed, Mag, Pha, PhaSy, old_Pha
						, old_PhaSy, PhaseDiff, dphi, freq, Y_out, size, Out,OutLen,num_win);
					///////////////////////

					random = rand() % (300 + 1);
					lshifter.setPeriod(from - random * step);

					printf("%d - %f\n", random, from - random * step);
					sprintf(index2, "%d", n);
					size_t len_path_tmp2 = strlen(path2) + strlen(label) + 1;
					char *temp2 = (char*)malloc(sizeof(char) * len_path_tmp);
					strcpy(temp2, path2);
					strcat(temp2, label);
					strcat(temp2, "_");
					size_t len = len_path_tmp + strlen(index) + 5;
					path_file2 = (char *)malloc(len * sizeof(char));
					strcpy(path_file2, temp2);
					strcat(path_file2, index2);
					strcat(path_file2, default_ext);
					//printf("path_F : %s \n", path_file2);
					/*if (!check_path(path_file2)) {
						continue;
					}*/
					fout = fopen(path_file2, "w");
					fprintf(fout, "%d\n", OutLen);
					for (int i = 0; i < OutLen; i++) {
						lshifter.tick(Out[i], fout);
					}
					fclose(fout);
					free(Out);
				}
				bound += numData;
			}
			i++;
			dem++;
		}
		++label_cur;
		i = 0;
	}
	fprintf(finf, "%d", sum);
	fclose(finf);
	free(cx_in);
	free(cx_out);
	free(Y);
	free(cx_out_i);
	kiss_fft_free(cfg);
	kiss_fft_free(cfg_i);
	free(Mag); free(old_PhaSy); free(freq); free(Y_out);
	free(Pha); free(PhaSy);  free(old_Pha); free(PhaseDiff); free(dphi);
}

void time_stretch(StkFloat *signal, int winsize, int fftsize, int shop, int ahop, float *CENTERFREQ, kiss_fft_cpx *Y, kiss_fft_cpx * cx_out_i,
	kiss_fft_cfg cfg_i, kiss_fft_cpx *cx_in, kiss_fft_cpx *cx_out, kiss_fft_cfg cfg, float *framed, float *Mag, float *Pha,
	float *PhaSy, float *old_Pha, float *old_PhaSy, float *PhaseDiff, float *dphi, float *freq, float *Y_out,
	int size, StkFloat *Out, int OutLen, int num_win) {
	for (int i = 0; i < fftsize; i++) {
		CENTERFREQ[i] = i * TWOPI*ahop / fftsize;
	}

	float temp;
	int dem = 0;
	int first = 0;
	int PosIn = 0;
	int PosOut = 0;

	for (int win_count = 0; win_count < num_win; win_count++) {
		if (first == 0) {

			for (int i = 0; i < winsize; i++) {
				framed[i] = signal[i] * HammingWindow(i, winsize);
				cx_in[i].r = framed[i];
				cx_in[i].i = 0;
			}
			kiss_fft(cfg, cx_in, cx_out);


			for (int j = 0; j < fftsize; j++) {
				Mag[j] = magnitude(cx_out[j].r, cx_out[j].i);
				Pha[j] = atan2(cx_out[j].i, cx_out[j].r);
				PhaSy[j] = Pha[j];
			}
			first = 1;
		}
		else {
			for (int i = PosIn, j = 0; i < PosIn + winsize; i++, j++) {
				framed[j] = signal[i] * HammingWindow(j, winsize);
				cx_in[j].r = framed[j];
				cx_in[j].i = 0;
			}
			kiss_fft(cfg, cx_in, cx_out);

			for (int j = 0; j < fftsize; j++) {
				Mag[j] = magnitude(cx_out[j].r, cx_out[j].i);
				Pha[j] = atan2(cx_out[j].i, cx_out[j].r);
				PhaseDiff[j] = Pha[j] - old_Pha[j];
				PhaseDiff[j] -= CENTERFREQ[j];
				dphi[j] = PhaseDiff[j] - (TWOPI)* round(PhaseDiff[j] / (TWOPI));		//FUCKMACRO
				freq[j] = (CENTERFREQ[j] + dphi[j]) / ahop;
				PhaSy[j] = old_PhaSy[j] + shop * freq[j];
			}
		}
		for (int i = 0; i < fftsize; i++) {
			Y[i].i = Mag[i] * sin(PhaSy[i]);
			Y[i].r = Mag[i] * cos(PhaSy[i]);
		}
		kiss_fft(cfg_i, Y, cx_out_i);
		for (int i = 0; i < fftsize; i++) {
			Y_out[i] = cx_out_i[i].r*HammingWindow(i, winsize);
		}
		for (int i = PosOut, j = 0; i < PosOut + winsize; i++, j++) {
			Out[i] += Y_out[j];
		}
		memcpy(old_Pha, Pha, sizeof(float) *fftsize);
		memcpy(old_PhaSy, PhaSy, sizeof(float) *fftsize);
		PosIn = PosIn + ahop;
		PosOut = PosOut + shop;
	}
}