/*
3次元_FDTD法による電磁界解析 ver. 2.01
From September, 2000;
Designed by Atsushi SAKAI;
supported by
Hiroya DAN (3D-_FDTD),
Yoshitaka WATANABE (Periodic Boundary Condition in 3D-_FDTD),
Hiroshi YAMADA(2D-_FDTD, PLRC),
Tomoki YONEHANA (2D-_FDTD PML, Non-Linear),
Toshihiko BABA (Photonic crystal bend model: April, 2001),
Kosuke MORITO (3D_Symmetry_Condition : October, 2003).
Takashi KAWASAKI (PCCW: 2007)
Koichiro YOSHIDA (Observation Area: 2008)
Norihiro ISHIKURA (October, 2012)
*/

#define _FDTD 1		// FDTD計算			0 : モデル掃き出し(プリプロセッサでコンパイルを変更させる)
//										1 : 計算実行

#define _BAND_CALCULATION 0			// 計算の種類 バンド計算
#define _PROPAGATION_CALCULATION 1	// 計算の種類 伝搬計算

#define _CALCULATION_TYPE _PROPAGATION_CALCULATION	// 計算の種類

#define _EXITATION_FUNC 1	// 励振関数の種類		0 : Gaussian
//													1 : CW

#define _PROGRAM_TEST 1		// プログラムの動作テスト	0: TEST(最終計算ステップ，出力ファイルを短く)
//															1: 本番

#define _MODEL_ALL_EPSILON 0 	// XY断面をモデル全体の出力	0: なし
//																1: あり

#define _CRT_SECURE_NO_WARNINGS //	警告を発生させないようにする

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//baba lab
//#include <direct.h>

//kuramitsu lab
#include<sys/stat.h>
#include<sys/types.h>



#include <stdlib.h>
#include <string>
#include "mpi.h"
#include "grobal_function.cpp"
#include "parameter.h"
#include "module0.h"

//サブルーティン
void file_open(char*);
void file_close();
void parameter(char*);
void initialize_matrix();
void modeling();
void set_epsilon();
void source_func();
void observation_func();
void calc_efield();
void calc_hfield();
void absorpt_bound_condition();
void saving_electric_field();
void output_field(char*);
void output_field_write(char *);
void output_model();

int main(int argc, char **argv){

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION
	double s_time, e_time;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	char time[9];
	int tag_send = 0, tag_recv = 0;
	int right, left;
	int namelen;

	MPI_Status status;

	// MPIによる通信の開始
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &isize);
	MPI_Comm_rank (MPI_COMM_WORLD, &irank);
	MPI_Get_processor_name (processor_name, &namelen);

	if (isize != ISIZE){
		printf ("MPIで設定した計算機の台数(%d)がプログラム中の値と一致しません．\n終了します\n", ISIZE);
		return 0;
	}

	printf ("%d分割並列処理スタート\n", isize);
	printf ("Process %d on %s\n", irank, processor_name);

	// 隣の計算機の番号の指定
	left = irank - 1;
	if(irank == IRANK_MIN){
		left = MPI_PROC_NULL;
	}
	right = irank + 1;
	if(irank == IRANK_MAX){
		right = MPI_PROC_NULL;
	}

	// dir_name (励振波長) の配列長だけ繰り返し
	for(int dir_count = 0; dir_count < (sizeof(dir_name) / sizeof(dir_name[0]) ); dir_count++){

		initialize_matrix(); 						// 配列の初期化
		modeling(); 								// モデルの設定
		file_open(dir_name[dir_count]); 			// ファイルを開く
		parameter(dir_name[dir_count]); 			// パラメータの設定と出力


		// 計算開始時刻の出力
		if (irank == IRANK_MIN){
			//_strtime(time);
			//fprintf(fpparameter, "Start Time:\t %s\n", time);
			//s_time = MPI_Wtime();
		}

		// 電磁界計算
		for(n = 1 ; n <= Nmax; n++){

			// 時間ステップ数の表示
			if(n % Ncut == 0){
				//_strtime(time);
				//printf("n = %d, \t\t", n);
				//printf("time = %s\n", time);
			}

			// 励振関数の設定
			source_func();

#if _FDTD

			// 一度同期をとる(同期はノード間で速度にばらつきが生じる作業)
			MPI_Barrier (MPI_COMM_WORLD);

			// 電界の計算
			calc_efield();

			// 吸収境界条件による端面の計算
			absorpt_bound_condition();

			// 一度同期をとる
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Sendrecv( &Ex[1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_send,
				&Ex[xmax][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ey[1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_send,
				&Ey[xmax][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ez[1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_send,
				&Ez[xmax][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ex[xmax-1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_send,
				&Ex[0][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ey[xmax-1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_send,
				&Ey[0][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ez[xmax-1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_send,
				&Ez[0][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

			// 電界の保存
			saving_electric_field();

			// 磁界の計算
			calc_hfield();

			MPI_Sendrecv( &Hy[xmax-1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_send,
				&Hy[0][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Hz[xmax-1][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, right, tag_send,
				&Hz[0][0][0], (ymax+1)*(zmax_ff+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

			// フィールドの出力
			output_field (dir_name[dir_count]);

			// ポインティングパワー計算と出力
#if _EXITATION_FUNC
#else
#endif

			// 一度同期をとる
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			if(n == 1) {
				observation_func(); 	// 観測点の設定
				output_model(); 		// モデルの出力
				set_epsilon(); 			// 誘電率の割り当て
			}
		}

		if (irank == IRANK_MIN){
			//_strtime(time);
			//fprintf(fpparameter, "End Time:\t %s\n", time); 	/*計算終了時刻の出力*/
			//時刻の出力
			//e_time = MPI_Wtime();
			//printf ("\ntime = %f\n", e_time - s_time);
		}

		file_close(); 			// ファイルを閉じる
	}

	//MPI_Finalize(); 			// MPIを終了する
#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
}



// 出力用ファイルを開く
void file_open(char* dir_name_def){
	char dir_name[40];
	char name_xy[40], name_yz[40], name_xz[40];

	sprintf(name_xy, "/Model_xy_%d.txt", irank);
	sprintf(name_yz, "/Model_yz_%d.txt", irank);
	sprintf(name_xz, "/Model_xz_%d.txt", irank);

	//baba lab
	//_mkdir(strcpy(dir_name, dir_name_def)); 		// 振り分けできるかテスト

	//kuramitsu lab
	mkdir(strcpy(dir_name, dir_name_def), 0755); 		// 振り分けできるかテスト



	if (irank == IRANK_MIN){
		fpparameter = fopen (strcat(strcpy(dir_name, dir_name_def), "/Parameter.txt"), "w");
		allmodel_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_xy.txt"), "w");
		fpallepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/All_Epsilon_xy.txt"), "w");
		allmodel_xz = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_xz.txt"), "w");
	}

	model_xy = fopen (strcat(strcpy(dir_name, dir_name_def), name_xy), "w"); 		// 振り分けできるかテスト
	model_yz = fopen (strcat(strcpy(dir_name, dir_name_def), name_yz), "w");
	model_xz = fopen (strcat(strcpy(dir_name, dir_name_def), name_xz), "w");

	fpepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_xy.txt"), "w");
	fpepsilony = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz.txt"), "w");
	fpepsilonz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx.txt"), "w");
	fpepsilony2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz2.txt"), "w");
	fpepsilonz2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx2.txt"), "w");

}


/*出力用ファイルを閉じる*/
void file_close(){

	if (irank == IRANK_MIN){
		fclose(fpparameter);
		fclose(fpallepsilonx);
	}

	fclose(model_xy);
	fclose(model_yz);
	fclose(model_xz);

	fclose(fpepsilonx);
	fclose(fpepsilony);
	fclose(fpepsilonz);
	fclose(fpepsilony2);
	fclose(fpepsilonz2);
}


// 計算用パラメータの設定と出力
void parameter(char* dir_name){

	if (irank == IRANK_MIN){

		fprintf(fpparameter, "XMAX_ALL = %d\n", XMAX_ALL);
		fprintf(fpparameter, "YMAX_ALL = %d\n", YMAX_ALL);
		fprintf(fpparameter, "ZMAX_ALL = %d\n", ZMAX_ALL);
		fprintf(fpparameter, "Nodes = %d\n", NODE);
		fprintf(fpparameter, "Cell Size [nm] = %d\n", CELL_SIZE);
		fprintf(fpparameter, "Time Step [s] = %e\n", dt);
		fprintf(fpparameter, "Final Time Step = %d\n", Nmax);
		fprintf(fpparameter, "Final Time [s] = %e\n", (double) Nmax * dt);
		fprintf(fpparameter, "\n");

		fprintf(fpparameter, "Upper Clad Index = %lf\n", n_clad);
		fprintf(fpparameter, "Slab Index = %lf\n", n_core);
		fprintf(fpparameter, "Upper Height [nm] = %d\n", CLAD_HEIGHT1);
		fprintf(fpparameter, "Slab Height [nm] = %d\n", SLAB_HEIGHT);
		fprintf(fpparameter, "Hole Pitch [nm] = %d\n", PITCH);
		fprintf(fpparameter, "Pitch Shift Max [nm] = %d\n", PITCH_SHIFT_MAX);
		fprintf(fpparameter, "Hole Diameter [nm] = %d\n", RADIUS*2);
		fprintf(fpparameter, "%d-row Shift [nm] = %d\n", 1, SX1);
		fprintf(fpparameter, "%d-row Shift [nm] = %d\n", 3, SX3);
		fprintf(fpparameter, "%Y-direction All Shift [nm] = %d\n", SY);
		fprintf(fpparameter, "Normal PCW Period = %d\n", NORM_PCW_PER);
		fprintf(fpparameter, "Chirped LSPCW Period = %d\n", CHIRP_3RD_LS_PER);
		fprintf(fpparameter, "Pitch Shift PCW Period = %d\n", PITCH_SHIFT_PER);
		fprintf(fpparameter, "Pitch Shift Chirp PCW Period = %d\n", PITCH_SHIFT_CHIRP_PER);
		fprintf(fpparameter, "LSPCW Period = %d\n", LSPCW_PER);
		fprintf(fpparameter, "Hole Column[y] = %d \n", intPcwWid);
		fprintf(fpparameter, "Hole Row[x] = %d \n", intPcwPer);
		fprintf(fpparameter, "Hole Start Coordinate = (%d, %d)\n", intPcwStartX, intPcwStartY);
		fprintf(fpparameter, "\n");

		fprintf(fpparameter, "Exctation [nm] = %d\n", EXCT_LEN);
		fprintf(fpparameter, "Exctation to Observation [nm] = %d\n", EXCT_OBSE_LEN);
		fprintf(fpparameter, "Observation to Wire [nm] = %d\n", OBSE_WIRE_LEN);
		fprintf(fpparameter, "Output Wire to Termination [nm] = %d\n", WIRE_OUTPUT_LEN);
		fprintf(fpparameter, "Termination Length [nm] = %d\n", WIRE_OUTPUT_OFFSET);
		if(intWireWid < 0){
			fprintf(fpparameter, "Wire Width [nm] = %d\n", (intWireWid + 1) * CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "Wire Width [nm] = %d\n", intWireWid * CELL_SIZE);
		}
		fprintf(fpparameter, "PCW Slab Termination Length [nm] = %d\n", PCW_SiSLAB_TERMINATION_LEN);
		if(PCW_SiSLAB_OFFSET < 0){
			fprintf(fpparameter, "PCW Slab Offset [nm] = %d\n", PCW_SiSLAB_OFFSET + CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "PCW Slab Offset [nm] = %d\n", PCW_SiSLAB_OFFSET);
		}
		if(PCW_WIDTH_CHIRP < 0){
			fprintf(fpparameter, "PCW Width Chirp [nm] = %d\n", PCW_WIDTH_CHIRP + CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "PCW Width Chirp [nm] = %d\n", PCW_WIDTH_CHIRP);
		}

		if(LSPCW_SHIFT_DESCRETE == FALSE){
			fprintf(fpparameter, "LSPCW_SHIFT_DESCRETE = FALSE\n");
		}
		else{
			fprintf(fpparameter, "LSPCW_SHIFT_DESCRETE = TRUE\n");
		}


		fprintf(fpparameter, "\n");
	}

	// 励振関数定数の設定
	lambda = atof(dir_name) * 1e-9;
	omega0 = 2.0*PI*C0/lambda;
	sigma = omega0 * delta_omega;
}

/*配列の初期化*/
void initialize_matrix(){

	//各ノードの座標
	if(irank != IRANK_MAX){
		xmax = XMAX;
		ymax = YMAX;
		zmax = ZMAX;
		zmax_ff = ZMAX_FF;
	}

	//最後のノードだけのりしろ不要なのでx方向に1セル小さい
	if(irank == IRANK_MAX){
		xmax = XMAX - 1;
		ymax = YMAX;
		zmax = ZMAX;
		zmax_ff = ZMAX_FF;
	}

	// 解析空間の最大値
	xmax_all = XMAX_ALL;
	ymax_all = YMAX_ALL;
	zmax_all = ZMAX_ALL;

	// 解析空間の中心座標
	x_cen = xmax/2;
	y_cen = ymax/2;
	z_cen = zmax/2;

	//モデルの中心と解析空間の中心は１セル分ずれているので要注意
	x_model_cen = x_cen + 1;
	y_model_cen = y_cen + 1;

	int x, y, z;

	for (x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax_ff+1; z++){
				// 電界
				Ex[x][y][z] = 0.0;
				Ey[x][y][z] = 0.0;
				Ez[x][y][z] = 0.0;

				// 磁界
				Hx[x][y][z] = 0.0;
				Hy[x][y][z] = 0.0;
				Hz[x][y][z] = 0.0;
			}
		}
	}

	for (x = 0; x < xmax_all; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax+1; z++){
				// 誘電率
				ALL_epsilonx[x][y][z] = EPSILON0;
				ALL_epsilony[x][y][z] = EPSILON0;
				ALL_epsilonz[x][y][z] = EPSILON0;
			}
		}
	}

	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax_ff+1; z++){
				// 誘電率(分割する必要がないような．．．)
				epsilonx[x][y][z] = EPSILON0;
				epsilony[x][y][z] = EPSILON0;
				epsilonz[x][y][z] = EPSILON0;
			}
		}
	}

	for(x = 0; x < xmax_all; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax+1; z++){
				// セルの目印
				ALL_cell[x][y][z] = CLAD;
			}
		}
	}


	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax_ff; z++){
				// セルの目印
				cell[x][y][z] = 0;
			}
		}
	}


	/****************************** Murの吸収境界条件 ******************************/

	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax_ff+1; z++){
			Exn2y00[x][z] = 0.0;
			Exn1y00[x][z] = 0.0;
			Exn2y01[x][z] = 0.0;
			Exn1y01[x][z] = 0.0;
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			Exn2z00[x][y] = 0.0;
			Exn1z00[x][y] = 0.0;
			Exn2z01[x][y] = 0.0;
			Exn1z01[x][y] = 0.0;
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Eyn2z00[x][y] = 0.0;
			Eyn1z00[x][y] = 0.0;
			Eyn2z01[x][y] = 0.0;
			Eyn1z01[x][y] = 0.0;
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax_ff+1; z++){
			Eyn2x00[y][z] = 0.0;
			Eyn1x00[y][z] = 0.0;
			Eyn2x01[y][z] = 0.0;
			Eyn1x01[y][z] = 0.0;
			Eyn2xm1[y][z] = 0.0;
			Eyn1xm1[y][z] = 0.0;
			Eyn2xm0[y][z] = 0.0;
			Eyn1xm0[y][z] = 0.0;
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax_ff; z++){
			Ezn2y00[x][z] = 0.0;
			Ezn1y00[x][z] = 0.0;
			Ezn2y01[x][z] = 0.0;
			Ezn1y01[x][z] = 0.0;
		}
	}
	for(y = 0; y <= ymax; y++){
		for(z = 0; z <= zmax_ff-1; z++){
			Ezn2x00[y][z] = 0.0;
			Ezn1x00[y][z] = 0.0;
			Ezn2x01[y][z] = 0.0;
			Ezn1x01[y][z] = 0.0;
			Ezn2xm1[y][z] = 0.0;
			Ezn1xm1[y][z] = 0.0;
			Ezn2xm0[y][z] = 0.0;
			Ezn1xm0[y][z] = 0.0;
		}
	}

	/****************************** Murの吸収境界条件 ******************************/
}


// モデルの設定
void modeling(){

	int n_temp; 		//屈折率の値保存用
	double epsilon_temp; 		//誘電率の値保存用

	/****************************** スラブの形成 ******************************/

	for(x = 0; x < xmax_all+1; x++){
		for(y = 0; y < ymax_all; y++){
			for(z = 0; z < zmax_all; z++){
				n_temp = CLAD;
				epsilon_temp = epsilon2;

				if(z < air_hc){			//空気層に設定
				}
				if(z >= air_hc && z < (air_hc + intCladHeight1) ){			//上部クラッドに設定
				}
				if(z >= (air_hc + intCladHeight1) && z < (air_hc + intCladHeight1 + intSlabHeigPer)){	// スラブに設定
					n_temp = CORE;
					epsilon_temp = epsilon1;
				}

				ALL_cell[x][y][z] = n_temp;
				ALL_epsilonx[x][y][z] = epsilon_temp;
				ALL_epsilony[x][y][z] = epsilon_temp;
				ALL_epsilonz[x][y][z] = epsilon_temp;
			}
		}
	}
	/****************************** スラブの形成 ******************************/


	/****************************** フォトニック結晶 ******************************/

	int s_x3; 		// チャープLSPCWのシフト量
	int s_x2;
	int s_x4;
	int z_end; 				// 円孔の終了座標


	if (intPcwPer == 0){
		intWirePer2 = intWireLen1 - 1;									// 出射COREスラブの開始点
		intWirePer3 = intWirePer2 + intWireLen2;						// 出射COREスラブの終了点
	}

	else{
		struct PNUM Pnum[100][10]; 	// 円柱の中心座標
		struct PNUM Pnum_Init[1][10]; 		// 円柱の標準格子定数による中心座標

		z_end = zmax_all; 		// 円孔が貫通している場合を考える
		Pnum_Init[0][0].Y = intPcwStartY; 	// 円柱を配置する最初のY座標

		if(intPcwWid % 2 == 1){		// if y:even
			Pnum_Init[0][intPcwWid-1].X = intWireLen1 + intPcwStartX + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
		}
		else{						// if y:odd 0.5Aずらす
			Pnum_Init[0][intPcwWid-1].X = intWireLen1 + intPcwStartX - 1;	// 配列の引数に使用するので-1
		}

		if(y != 0){
			Pnum_Init[0][intPcwWid-1].Y = Pnum_Init[0][0].Y + intPitchY * (intPcwWid - 1); 		//ｙ座標を(root3)/2*intPitchXだけずらす
		}

		for(z = 0; z < z_end; z++){
			for(y = 0; y < intPcwWid; y++){

				int input_NormPcw_Xend;
				int input_PitchShiftPcw_Xend;
				int input_PitchShiftChirpPcw_Xend;
				int input_Chirp_Ls_Xend;
				int Lspcw_Xend;
				int output_Chirp_Ls_Xend;
				int output_PitchShiftChirpPcw_Xend;
				int output_PitchShiftPcw_Xend;
				int output_PCW_Xend;

				if (intNormPcwPer > 0){
					input_NormPcw_Xend = intNormPcwPer;
				}
				else{
					input_NormPcw_Xend = 0;
				}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				input_PitchShiftPcw_Xend = input_NormPcw_Xend + intPitchShiftPcwPer + intPitchShiftChirpPcwPer;
				input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/

				if (intChirp3rdLsPer > 0){
					input_Chirp_Ls_Xend = input_PitchShiftChirpPcw_Xend + (intChirp3rdLsPer);
				}
				else{
					input_Chirp_Ls_Xend = input_PitchShiftChirpPcw_Xend;
				}
				Lspcw_Xend = input_Chirp_Ls_Xend + intLspcwPer + 1;
				if (intChirp3rdLsPer > 0){
					output_Chirp_Ls_Xend = Lspcw_Xend + (intChirp3rdLsPer - 1);
				}
				else{
					output_Chirp_Ls_Xend = Lspcw_Xend;
				}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend;
				output_PitchShiftPcw_Xend = output_PitchShiftChirpPcw_Xend + intPitchShiftPcwPer + intPitchShiftChirpPcwPer;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/

				output_PCW_Xend = output_Chirp_Ls_Xend + intNormPcwPer;

				int intPitchShiftX = INT_DIV (PITCH_SHIFT_MAX, CELL_SIZE);
				int intPitchShiftY = (INT_DIV((PITCH_SHIFT_MAX * sqrt(3.0)/2 + 0.5), CELL_SIZE));
				int intPCWwidthChirp = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
				int intPreviousPCWwidthOffset;
				int intNowPCWwidthOffset;

				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				double  dblPitchShiftChirpY;
				//double dblPitchShiftChirpX, dblPitchShiftChirpY, dblPitchShiftChirpY2;
				int intPitchShiftChirpX, intPitchShiftChirpY;
				//int intPitchShiftChirpX, intPitchShiftChirpX2, intPitchShiftChirpY;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/


				/****************************** LSPCW ******************************/
				for(x = 0; x < intPcwPer; x++){
					int y2, y_poo, y_poo2;

					s_x3 = 0;
					s_x2 = 0;
					s_x4 = 0;
					// 入射 通常PCW
					if (x < input_NormPcw_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][y2].Y;
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;
							}
						}
					}

					// 入射 格子定数変化PCW
					else if (x < input_PitchShiftPcw_Xend && x >= input_NormPcw_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchShiftY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchShiftX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchShiftX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路1列目シフト構造(2013/7/12) ********************/
								if (y2 != intPcwWid - 1){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								/******************** 導波路1列目シフト構造(2013/7/12) ********************/

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							/******************** 格子シフト量チャープ(2013/7/19) ********************/
							intPitchShiftChirpX = intPitchShiftX;
							intPitchShiftChirpY = 0;

							if (x >= intPitchShiftPcwPer){
								intPitchShiftChirpX -= (int) ((intPitchShiftX - intPitchX) * ((x - intPitchShiftPcwPer - 1) / (double) (intPitchShiftChirpPcwPer - 1)));

								if (x > intPitchShiftPcwPer){
									dblPitchShiftChirpY = (intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
								}
							}

							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								if(x < intPitchShiftPcwPer + intPitchShiftChirpPcwPer - 1){
									if (y2 < intPcwWid-1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchShiftY - dblPitchShiftChirpY * (x - intPitchShiftPcwPer)) * (intPcwWid-1 - y2));
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][y2].Y;
									}
								}
								else{
									if (y2 < intPcwWid-1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][y2].Y;
									}
								}
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
							}
							/******************** 格子シフト量チャープ(2013/7/19) ********************/

						}
						if (LSPCW_SHIFT_DESCRETE == FALSE){
							// 2列目格子シフト
							if (y == intPcwWid - 2){
								s_x2 = INT_DIV(SX2, CELL_SIZE);
							}
							// 3列目格子シフト
							if (y == intPcwWid - 3){
								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								s_x3 = INT_DIV(SX3, CELL_SIZE);
								/******************** 格子シフト量チャープ(2013/7/19) ********************/

							}
							if (y == intPcwWid - 4){
								s_x4 = INT_DIV(SX4, CELL_SIZE);
							}
						}
					}

					// 入射 チャープLSPCW
					else if (x >= input_NormPcw_Xend && x < input_Chirp_Ls_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][y2].Y;
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;

							}

							if (y == intPcwWid - 3){
								if (intChirp3rdLsPer == 0){
									s_x3 = 0;
								}
								else{
									s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (x - input_NormPcw_Xend);
								}
							}
							if (y == intPcwWid - 2){
								if (intChirp3rdLsPer == 0){
									s_x2 = 0;
								}
								else{
									s_x2 = INT_DIV (intSx2Per, intChirp3rdLsPer) * (x - input_NormPcw_Xend);
								}
							}
							if (y == intPcwWid - 4){
								if (intChirp3rdLsPer == 0){
									s_x4 = 0;
								}
								else{
									s_x4 = INT_DIV (intSx4Per, intChirp3rdLsPer) * (x - input_NormPcw_Xend);
								}
							}
						}
					}


					// 出射 格子定数変化PCW
					else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){

						// 格子定数変化PCWとの入射接続部
						if (x == output_PitchShiftChirpPcw_Xend){
							y_poo = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){

								/******************** 導波路1列目シフト構造(2013/7/12) ********************/
								if (y2 != intPcwWid - 1){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								/******************** 導波路1列目シフト構造(2013/7/12) ********************/

								y_poo++;
							}
						}
						/******************** 格子シフト量チャープ(2013/7/19) ********************/
						intPitchShiftChirpX = intPitchX;
						intPitchShiftChirpY = 0;
						dblPitchShiftChirpY = 0;

						if (x >= output_PitchShiftChirpPcw_Xend){
							intPitchShiftChirpX += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend) / (double) (intPitchShiftChirpPcwPer - 1)));
							dblPitchShiftChirpY = (intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
						}

						for (y2 = intPcwWid-1; y2 >= 0; y2--){
							if (x != output_PitchShiftChirpPcw_Xend){
								if (y2 < intPcwWid-1){
									if (y2 % 2 == 1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend)) * (intPcwWid-1 - y2));
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend + 1)) * (intPcwWid-1 - y2));
									}
								}
								else{
									Pnum[x][y2].Y = Pnum[x-1][y2].Y;
								}
							}
							else{
								if (y2 < intPcwWid-1){
									if (y2 % 2 == 1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);
									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend + 1)) * (intPcwWid-1 - y2));
									}
								}
								else{
									Pnum[x][y2].Y = Pnum[x-1][y2].Y;
								}
							}

							Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						}

						/******************** 格子シフト量チャープ(2013/7/19) ********************/

						if (LSPCW_SHIFT_DESCRETE == FALSE){
							// 3列目格子シフト
							if (y == intPcwWid - 2){
								s_x2 = INT_DIV(SX2, CELL_SIZE);
							}
							if (y == intPcwWid - 4){
								s_x4 = INT_DIV(SX4, CELL_SIZE);
							}
							if (y == intPcwWid - 3){

								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								s_x3 = INT_DIV(SX3, CELL_SIZE);
								/******************** 格子シフト量チャープ(2013/7/19) ********************/
							}
						}
					}

					// 出射 チャープLSPCW
					else if (x >= Lspcw_Xend && x < output_Chirp_Ls_Xend){
						for (y2 = intPcwWid-1; y2 >= 0; y2--){
							Pnum[x][y2].Y = Pnum[x-1][y2].Y;
							Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;

						}
						if (y == intPcwWid - 3){
							if (intChirp3rdLsPer == 0){
								s_x3 = 0;
							}
							else{
								s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x);
							}
						}
						if (y == intPcwWid - 4){
							if (intChirp3rdLsPer == 0){
								s_x4 = 0;
							}
							else{
								s_x4 = INT_DIV (intSx4Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x);
							}
						}
						if (y == intPcwWid - 2){
							if (intChirp3rdLsPer == 0){
								s_x2 = 0;
							}
							else{
								s_x2 = INT_DIV (intSx2Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x);
							}
						}
					}


					// 出射 通常PCW
					else if (x >= output_Chirp_Ls_Xend && x < output_PCW_Xend){
						if (x == 0){
							flag_2r = 1;
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								if(y2 % 2 == 1){		// if y:even
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][y2].Y;
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;
							}
						}
					}

					// 通常格子定数PCW or LSPCW
					else{
						if (x != 0){
							y_poo = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (intPitchY * y_poo);
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;
								y_poo++;
							}
						}

						// 3列目格子シフト
						if (y == intPcwWid - 3){
							s_x3 = intSx3Per;
						}
						if (y == intPcwWid - 2){
							s_x2 = intSx2Per;
						}
						if (y == intPcwWid - 4){
							s_x4 = intSx4Per;
						}
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// 導波路幅チャープ
					if (PCW_WIDTH_CHIRP != 0){
						if (x < input_PitchShiftPcw_Xend){
							////////入射シフト量チャープここから(田村)
							if(intChirp2ndLsPer != 0){
							//3列目
							if (y == intPcwWid - 3){
								if (intChirp2ndLsPer == 0){
									s_x3 = 0;
								}
								else if (INT_DIV (intSx3Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend)<INT_DIV(SX3,CELL_SIZE)){
									s_x3 = INT_DIV (intSx3Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend);
								}
								else{
									s_x3 = INT_DIV(SX3,CELL_SIZE);
								}
							}
							//2列目
							if (y == intPcwWid - 2){
								if (intChirp2ndLsPer == 0){
									s_x2 = 0;
								}
								else if (INT_DIV (intSx2Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend)<INT_DIV(SX2,CELL_SIZE)){
									s_x2 = INT_DIV (intSx2Per, intChirp2ndLsPer) * (x +10 - input_PitchShiftPcw_Xend);
								}
								else{
									s_x2 = INT_DIV(SX2,CELL_SIZE);
								}
							}
							}
							////////入射シフト量チャープここまで
							if (x == 0){
								intNowPCWwidthOffset = intPCWwidthChirp;
								for (y2 = intPcwWid-1; y2 >= 0; y2--){
									Pnum[x][y2].Y -= intNowPCWwidthOffset;
								}
								intPreviousPCWwidthOffset = intNowPCWwidthOffset;
							}
							else{
								// PCW_WIDTH_CHIRP用微調整
								intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * (x / (double) (input_PitchShiftPcw_Xend-1)));
								//intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * ((x + 1) / (double) (input_PitchShiftPcw_Xend-1)));
								if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
									}
									intPreviousPCWwidthOffset = intNowPCWwidthOffset;
								}
							}
						}
						else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
							////////出射シフト量チャープここから(田村)
							if(intChirp2ndLsPer != 0){
							//3列目
							if (y == intPcwWid - 3){
								if (intChirp2ndLsPer == 0){
									s_x3 = 0;
								}
								else if( INT_DIV (intSx3Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2)<INT_DIV(SX3,CELL_SIZE)){
									s_x3 = INT_DIV (intSx3Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2);
								}
								else{
									s_x3 = INT_DIV(SX3,CELL_SIZE);
								}
							}
							//2列目
							if (y == intPcwWid - 2){
								if (intChirp2ndLsPer == 0){
									s_x2 = 0;
								}
								else if( INT_DIV (intSx2Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2)<INT_DIV(SX2,CELL_SIZE)){
									s_x2 = INT_DIV (intSx2Per, intChirp2ndLsPer) * (output_PitchShiftPcw_Xend - x-2);
								}
								else{
									s_x2 = INT_DIV(SX2,CELL_SIZE);
								}
							}
							}
						//////////出射シフト量チャープここまで

							if (x == output_PitchShiftChirpPcw_Xend){
								intPreviousPCWwidthOffset = 0;

								intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
								if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										if (y2 % 2 == 0){
											Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
										}
									}
									intPreviousPCWwidthOffset = intNowPCWwidthOffset;
								}
								intPreviousPCWwidthOffset = 0;
							}
							else{
								// PCW_WIDTH_CHIRP用微調整
								for (y2 = intPcwWid-1; y2 >= 0; y2--){
									if (y2 % 2 == 1){
										intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
									}
									else{
										intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
									}
									if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
										Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
									}
								}

								intPreviousPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)) + 0.9);
							}
						}
					}
					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


					/******************** 格子シフト量チャープ(2013/7/19) ********************/
					if (PITCH_SHIFT_MAX > PITCH && intNormPcwPer == 0 && intPitchShiftPcwPer + intPitchShiftChirpPcwPer != 0){
					}
					/******************** 格子シフト量チャープ(2013/7/19) ********************/

					flag_2r = 0;
					if(x < PITCH_SHIFT_PER && x > PITCH_SHIFT_PER+LSPCW_PER){
						flag_2r = 0;
					}else if(x%2 == 0 && x >= PITCH_SHIFT_PER && x <= PITCH_SHIFT_PER+LSPCW_PER){
						flag_2r = 1;
					}else if(x%2 == 1 && x >= PITCH_SHIFT_PER && x <= PITCH_SHIFT_PER+LSPCW_PER){
						flag_2r = 2;
					}

					if (SX2 == 0 && SX4 == 0)
					mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
					if (SX3 == 0 && SX4 == 0)
					mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
					if (SX2 == 0 && SX3 == 0)
					mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);


					// 出射CORE細線導波路との接続部分の調整
					if(x == intPcwPer - 1){
						// 格子定数変化PCWあり
						/******************** 格子シフト量チャープ(2013/7/19) ********************/
						if (intNormPcwPer == 0 && intPitchShiftPcwPer + intPitchShiftChirpPcwPer != 0){
							/******************** 格子シフト量チャープ(2013/7/19) ********************/

							if(y % 2 == 1){
								if (Pnum[x][y].X > 0){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									if (PITCH_SHIFT_MAX == PITCH){
										//入ってる
										Pnum[x][y].X += intPitchX;
									}
									else{
										Pnum[x][y].X += intPitchShiftX; 		// intPitchXだけ+X座標に置く
									}
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								else if (Pnum[x-1][y].X > 0){
									Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-2][y].X > 0){
									Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}

								int poo = 0;

								if (PCW_WIDTH_CHIRP != 0){
									//入ってる
									// PCW_WIDTH_CHIRP用微調整
									poo = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
									Pnum[x][y].Y -= INT_DIV(poo, (intPitchShiftPcwPer-1));
								}


								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								if (PITCH_SHIFT_MAX != PITCH){
									Pnum[x][y].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchShiftY * (intPcwWid-1 - y) - intPCWwidthChirp;
								}
								/******************** 格子シフト量チャープ(2013/7/19) ********************/

								if (SX2 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
								if (SX3 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
								if (SX2 == 0 && SX3 == 0)
									mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);
							}
						}
						else{
							if(y % 2 == 1){

								if (Pnum[x][y].X > 0){
									Pnum[x][y].X += intPitchX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-1][y].X > 0){
									Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-2][y].X > 0){
									Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}

								if (PCW_WIDTH_CHIRP != 0){
									// PCW_WIDTH_CHIRP用微調整
									int poo;
									poo = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
									Pnum[x][y].Y -= INT_DIV(poo, (intPitchShiftPcwPer-1));
								}

								if (SX2 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
								if (SX3 == 0 && SX4 == 0)
									mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
								if (SX2 == 0 && SX3 == 0)
									mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);
							}
						}
					}
				}
				/****************************** LSPCW ******************************/

			}
		}
		// CORE細線導波路の位置を計算
		intWirePer2 = Pnum[intPcwPer-1][intPcwWid-1].X + intRadius;		// 出射COREスラブの開始点
		intWirePer3 = intWirePer2 + intWireLen2;						// 出射COREスラブの終了点

	}
	/****************************** フォトニック結晶 ******************************/

	/****************************** 入出射細線導波路 ******************************/
	int intPcwSislabOffset;

	// 全面スラブになっているので，細線以外の部分を空気に変更
	if (PCW_SiSLAB_OFFSET != 0){
		intPcwSislabOffset = INT_DIV(PCW_SiSLAB_OFFSET, CELL_SIZE);
	}
	else{
		intPcwSislabOffset = 0;
	}

	for (z = zmax_all - intSlabHeigPer; z < (zmax_all + 1); z++){
		for (y = 0; y < ymax_all - intWireWid_2; y++){

			// 入射
			if (PCW_SiSLAB_OFFSET != 0){
			}
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
				ALL_cell[x][y][z] = CLAD;
				ALL_epsilonx[x][y][z] = epsilon2;
				ALL_epsilony[x][y][z] = epsilon2;
				ALL_epsilonz[x][y][z] = epsilon2;
			}

			// 出射
			if (PCW_SiSLAB_OFFSET != 0){
			}
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++){
				ALL_cell[x][y][z] = CLAD;
				ALL_epsilonx[x][y][z] = epsilon2;
				ALL_epsilony[x][y][z] = epsilon2;
				ALL_epsilonz[x][y][z] = epsilon2;
			}
		}

	}
	/****************************** 入出射細線導波路 ******************************/



	/****************************** 対称境界部分の誘電率の設定 ******************************/

	for(x = 0; x < xmax_all+1; x++){
		for(z = 0; z < zmax_all+1; z++){
			ALL_cell[x][ymax][z] = ALL_cell[x][ymax-1][z];
			ALL_epsilonx[x][ymax][z] = ALL_epsilonx[x][ymax-1][z];
			ALL_epsilony[x][ymax][z] = ALL_epsilony[x][ymax-1][z];
			ALL_epsilonz[x][ymax][z] = ALL_epsilonz[x][ymax-1][z];
		}
	}

	for(x = 0; x < xmax_all+1; x++){
		for(y = 0; y < ymax_all+1; y++){
			ALL_cell[x][y][zmax] = ALL_cell[x][y][zmax-1];
			ALL_epsilonx[x][y][zmax] = ALL_epsilonx[x][y][zmax-1];
			ALL_epsilony[x][y][zmax] = ALL_epsilony[x][y][zmax-1];
			ALL_epsilonz[x][y][zmax] = ALL_epsilonz[x][y][zmax-1];
		}
	}

	/****************************** 対称境界部分の誘電率の設定 ******************************/



	/****************************** 各ノードにモデルを分割 ******************************/

	if(NODE < 3){
		if(irank != IRANK_MAX){
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < Z_AIR; z++){
						cell[x][y][z] = AIR;
						epsilonx[x][y][z] = epsilon0;
						epsilony[x][y][z] = epsilon0;
						epsilonz[x][y][z] = epsilon0;
					}
					for(z = Z_AIR; z < zmax_ff+1; z++){
						cell[x][y][z] = ALL_cell[irank*(xmax)+x][y][z-Z_AIR];
						epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax-1)+x][y][z-Z_AIR];
						epsilony[x][y][z] = ALL_epsilony[irank*(xmax-1)+x][y][z-Z_AIR];
						epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax-1)+x][y][z-Z_AIR];
					}
				}
			}

		}
		else{
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < Z_AIR; z++){
						cell[x][y][z] = AIR;
						epsilonx[x][y][z] = epsilon0;
						epsilony[x][y][z] = epsilon0;
						epsilonz[x][y][z] = epsilon0;
					}
					for(z = Z_AIR; z < zmax_ff+1; z++){
						cell[x][y][z] = ALL_cell[irank*(xmax)+x][y][z-Z_AIR];
						epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax)+x][y][z-Z_AIR];
						epsilony[x][y][z] = ALL_epsilony[irank*(xmax)+x][y][z-Z_AIR];
						epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax)+x][y][z-Z_AIR];
					}
				}
			}
		}
	}else{
		if(irank == IRANK_MAX){
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < Z_AIR; z++){
						cell[x][y][z] = AIR;
						epsilonx[x][y][z] = epsilon0;
						epsilony[x][y][z] = epsilon0;
						epsilonz[x][y][z] = epsilon0;
					}
					for(z = Z_AIR; z < zmax_ff+1; z++){
						cell[x][y][z] = ALL_cell[2*(xmax)+x][y][z-Z_AIR];
						epsilonx[x][y][z] = ALL_epsilonx[2*(xmax)+x][y][z-Z_AIR];
						epsilony[x][y][z] = ALL_epsilony[2*(xmax)+x][y][z-Z_AIR];
						epsilonz[x][y][z] = ALL_epsilonz[2*(xmax)+x][y][z-Z_AIR];
					}
				}
			}
		}else if(irank == IRANK_MIN){
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < Z_AIR; z++){
						cell[x][y][z] = AIR;
						epsilonx[x][y][z] = epsilon0;
						epsilony[x][y][z] = epsilon0;
						epsilonz[x][y][z] = epsilon0;
					}
					for(z = Z_AIR; z < zmax_ff+1; z++){
						cell[x][y][z] = ALL_cell[x][y][z-Z_AIR];
						epsilonx[x][y][z] = ALL_epsilonx[x][y][z-Z_AIR];
						epsilony[x][y][z] = ALL_epsilony[x][y][z-Z_AIR];
						epsilonz[x][y][z] = ALL_epsilonz[x][y][z-Z_AIR];
					}
				}
			}
		}else{
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < Z_AIR; z++){
						cell[x][y][z] = AIR;
						epsilonx[x][y][z] = epsilon0;
						epsilony[x][y][z] = epsilon0;
						epsilonz[x][y][z] = epsilon0;
					}
					for(z = Z_AIR; z < zmax_ff+1; z++){
						cell[x][y][z] = ALL_cell[(xmax-1)+x][y][z-Z_AIR];
						epsilonx[x][y][z] = ALL_epsilonx[(xmax-1)+x][y][z-Z_AIR];
						epsilony[x][y][z] = ALL_epsilony[(xmax-1)+x][y][z-Z_AIR];
						epsilonz[x][y][z] = ALL_epsilonz[(xmax-1)+x][y][z-Z_AIR];
					}
				}
			}

		}


	}



	/*if(NODE < 3){
		if(irank != IRANK_MAX){
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < zmax+1; z++){
						cell[x][y][z] = ALL_cell[irank*(xmax-1)+x][y][z];
						epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax-1)+x][y][z];
						epsilony[x][y][z] = ALL_epsilony[irank*(xmax-1)+x][y][z];
						epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax-1)+x][y][z];
					}
				}
			}
		}
		else{
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < zmax+1; z++){
						cell[x][y][z] = ALL_cell[irank*(xmax)+x][y][z];
						epsilonx[x][y][z] = ALL_epsilonx[irank*(xmax)+x][y][z];
						epsilony[x][y][z] = ALL_epsilony[irank*(xmax)+x][y][z];
						epsilonz[x][y][z] = ALL_epsilonz[irank*(xmax)+x][y][z];
					}
				}
			}
		}
	}else{
		if(irank == IRANK_MIN){
			if(NUM_EX < 3){
				for(x = 0; x < xmax+1; x++){
					for(y = 0; y < ymax+1; y++){
						for(z = 0; z < zmax+1; z++){
							cell[x][y][z] = ALL_cell[x][y][z];
							epsilonx[x][y][z] = ALL_epsilonx[x][y][z];
							epsilony[x][y][z] = ALL_epsilony[x][y][z];
							epsilonz[x][y][z] = ALL_epsilonz[x][y][z];
						}
					}
				}
			}else{
				for(x = 0; x < X_LSPCW-1; x++){
					for(y = 0; y < ymax+1; y++){
						for(z = 0; z < zmax+1; z++){
							cell[x][y][z] = ALL_cell[x][y][z];
							epsilonx[x][y][z] = ALL_epsilonx[x][y][z];
							epsilony[x][y][z] = ALL_epsilony[x][y][z];
							epsilonz[x][y][z] = ALL_epsilonz[x][y][z];
						}
					}
				}
				for(x = X_LSPCW-1; x < xmax+1; x++){
					for(y = 0; y < ymax+1; y++){
						for(z = 0; z < zmax+1; z++){
							cell[x][y][z] = ALL_cell[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
							epsilonx[x][y][z] = ALL_epsilonx[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
							epsilony[x][y][z] = ALL_epsilony[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
							epsilonz[x][y][z] = ALL_epsilonz[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
						}
					}
				}
			}

		}else if(irank == IRANK_MAX){
			if(NUM_EX < 3){
				for(x = 0; x < xmax+1; x++){
					for(y = 0; y < ymax+1; y++){
						for(z = 0; z < zmax+1; z++){
							cell[xmax-x][y][z] = ALL_cell[xmax_all-x][y][z];
							epsilonx[xmax-x][y][z] = ALL_epsilonx[xmax_all-x][y][z];
							epsilony[xmax-x][y][z] = ALL_epsilony[xmax_all-x][y][z];
							epsilonz[xmax-x][y][z] = ALL_epsilonz[xmax_all-x][y][z];
						}
					}
				}
			}else{
				for(x = 0; x < (X_LSPCW-1)*(NUM_EX-1)+1; x++){
					for(y = 0; y < ymax+1; y++){
						for(z = 0; z < zmax+1; z++){
							cell[x][y][z] = ALL_cell[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
							epsilonx[x][y][z] = ALL_epsilonx[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
							epsilony[x][y][z] = ALL_epsilony[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
							epsilonz[x][y][z] = ALL_epsilonz[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
						}
					}
				}
				for(x = (X_LSPCW-1)*(NUM_EX-1)+1; x < xmax+1; x++){
					for(y = 0; y < ymax+1; y++){
						for(z = 0; z < zmax+1; z++){
							cell[x][y][z] = ALL_cell[2*(X_LSPCW-1)+(x-(X_LSPCW-1)*(NUM_EX-1))][y][z];
							epsilonx[x][y][z] = ALL_epsilonx[2*(X_LSPCW-1)+(x-(X_LSPCW-1)*(NUM_EX-1))][y][z];
							epsilony[x][y][z] = ALL_epsilony[2*(X_LSPCW-1)+(x-(X_LSPCW-1)*(NUM_EX-1))][y][z];
							epsilonz[x][y][z] = ALL_epsilonz[2*(X_LSPCW-1)+(x-(X_LSPCW-1)*(NUM_EX-1))][y][z];
						}
					}
				}

			}
		}else{
			for(x = 0; x < xmax+1; x++){
				for(y = 0; y < ymax+1; y++){
					for(z = 0; z < zmax+1; z++){
						cell[x][y][z] = ALL_cell[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
						epsilonx[x][y][z] = ALL_epsilonx[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
						epsilony[x][y][z] = ALL_epsilony[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
						epsilonz[x][y][z] = ALL_epsilonz[x%(X_LSPCW-1)+(X_LSPCW-1)][y][z];
					}
				}
			}

		}

	}*/





	/****************************** 各ノードにモデルを分割 ******************************/


	/****************************** 共通パラメータの設定 ******************************/

	// 励振点，観測面の設定 (XMAXは "のりしろ" 部分を含めていることに注意)
	intExctPortNum = intExctLen / (XMAX - 1);
	intObseInPortNum = intObseLen1 / (XMAX - 1);
	intObseOutPortNum = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) / (XMAX - 1);
	if (NODE % 2 != 0){
		intObseCenPortNum = XMAX_ALL / 2 / (XMAX - 1);	// 奇数個並列計算のとき
	}
	else{
		intObseCenPortNum = XMAX_ALL / 2 / (XMAX - 1) - 1;	// 偶数個並列計算のとき
	}

	intExctLenPart = intExctLen % (XMAX - 1) - 1;		// 配列の引数に使用するので-1
	intObseLenPart1 = intObseLen1 % (XMAX - 1) - INT_DIV(intObseInter, 2) - 1;
	intObseLenPart2 = intObseLenPart1 + intObseInter;
	intObseLenPart3 = intObseLen1 % (XMAX - 1);
	intObseLenPart4 = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1) - INT_DIV(intObseInter, 2);
	intObseLenPart5 = intObseLenPart4 + intObseInter;
	intObseLenPart6 = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1);
	if (NODE % 2 != 0){
		intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1) - 10;		// 奇数個並列計算のとき
	}
	else{
		intObseLenPart7 = (XMAX - 1) - 10;		// 偶数個並列計算のとき
	}
	/****************************** 共通パラメータの設定 ******************************/
}



//誘電率の割り当て
void set_epsilon(){

	//誘電率分布の出力(モデルの確認)
	int tag1 = 1;

#if _FDTD

	/****************************** 計算実行時 ******************************/
	int node;

	MPI_Status status;

	//XY平面
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			epsilon_xy[x][y] = epsilonx[x][y][intSlabCen-1];
		}
	}
	if(irank == IRANK_MIN){
		for(x = 0; x<xmax; x++){
			for(y = 0; y < ymax+1; y++){
				fprintf (fpallepsilonx, "%e\t", epsilon_xy[x][y]);
				fprintf (fpepsilonx, "%e\t", epsilon_xy[x][y]);
			}
			fprintf (fpallepsilonx, "\n");
			fprintf (fpepsilonx, "\n");
		}
	}
	if(irank != IRANK_MIN){

		MPI_Send (&epsilon_xy[0][0], (XMAX+1)*(YMAX+1), MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD);

		for(x = 1; x<xmax; x++){
			for(y = 0; y < ymax+1; y++){
				fprintf (fpepsilonx, "%e\t", epsilon_xy[x][y]);
			}
			fprintf (fpepsilonx, "\n");
		}
	}
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){

			MPI_Recv (&epsilon_xy[0][0], (XMAX+1)*(YMAX+1), MPI_DOUBLE, node, tag1, MPI_COMM_WORLD, &status);

			if(node == IRANK_MAX){
				for(x = 1; x < xmax-1; x++){
					for(y = 0; y < ymax+1; y++){
						fprintf(fpallepsilonx, "%e\t", epsilon_xy[x][y]);
					}
					fprintf(fpallepsilonx, "\n");
				}
			}
			else{
				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax+1; y++){
						fprintf (fpallepsilonx, "%e\t", epsilon_xy[x][y]);
					}
					fprintf (fpallepsilonx, "\n");
				}
			}
		}
	}
	/****************************** 計算実行時 ******************************/
#else

	/****************************** モデル確認時 ******************************/
	char fname[40],dir_name[50];	//ファイル名格納変数

	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			epsilon_xy[x][y] = epsilonx[x][y][intSlabCen-1];
		}
	}
	for(x = 0; x<xmax; x++){
		for(y = 0; y<ymax+1; y++){
			fprintf (fpallepsilonx, "%e\t", epsilon_xy[x][y]);
			fprintf (fpepsilonx, "%e\t", epsilon_xy[x][y]);
		}
		fprintf (fpallepsilonx, "\n");
		fprintf (fpepsilonx, "\n");
	}

#if _PROGRAM_TEST
#if	_MODEL_ALL_EPSILON
	//for(z = 0; z < zmax+1; z++){
	for(z = 0; z < zmax+1; z++){
		//for(z = intCladHeight1 - 1; z < zmax+1; z++){
		sprintf(fname,"/AllEpsilon_x_%d.txt",z);
		fpAllEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");
		sprintf(fname,"/Epsilon_x_%d.txt",z);
		fpEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");

		for(x=0; x < xmax; x++){
			for(y=0; y < ymax; y++){
				fprintf (fpAllEpsilon, "%e\t", epsilonx[x][y][z]);
				fprintf (fpEpsilon, "%e\t", epsilonx[x][y][z]);
			}
			fprintf (fpAllEpsilon, "\n");
			fprintf (fpEpsilon, "\n");
		}

		fclose(fpAllEpsilon);
		fclose(fpEpsilon);
	}
#endif
#else
	for(z = 0; z < zmax+1; z++){
		sprintf(fname,"/AllEpsilon_x_%d.txt",z);
		fpAllEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");
		sprintf(fname,"/Epsilon_x_%d.txt",z);
		fpEpsilon = fopen(strcat(strcpy(dir_name, "FOR_TEST"),fname),"w");

		for(x=0; x < xmax; x++){
			for(y=0; y < ymax+1; y++){
				fprintf (fpAllEpsilon, "%e\t", epsilonx[x][y][z]);
				fprintf (fpEpsilon, "%e\t", epsilonx[x][y][z]);
			}
			fprintf (fpAllEpsilon, "\n");
			fprintf (fpEpsilon, "\n");
		}

		fclose(fpAllEpsilon);
		fclose(fpEpsilon);
	}
#endif

	/****************************** モデル確認時 ******************************/
#endif

	//YZ平面
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax_ff+1; z++){
			epsilon_yz[y][z] = epsilony[intObseLenPart1][y][z];
		}
	}
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax_ff+1; z++){
			fprintf(fpepsilony, "%e\t", epsilon_yz[y][z]);
		}
		fprintf(fpepsilony, "\n");
	}

	//ZX平面 (Y:境界面)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff+1; z++){
			epsilon_zx[x][z] = epsilonz[x][ymax][z];
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff+1; z++){
			fprintf(fpepsilonz, "%e\t", epsilon_zx[x][z]);
		}
		fprintf(fpepsilonz, "\n");
	}
	//ZX平面 (Y:中心)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff+1; z++){
			epsilon_zx2[x][z] = epsilonz[x][ymax/2][z];
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff+1; z++){
			fprintf(fpepsilonz2, "%e\t", epsilon_zx2[x][z]);
		}
		fprintf(fpepsilonz2, "\n");
	}


	// ファイルポインタを閉じる
	if (irank == IRANK_MIN){
		fclose(fpallepsilonx);
	}
	fclose(fpepsilonx);
	fclose(fpepsilony);
	fclose(fpepsilonz);

}


// 励振関数
void source_func(){

	int x, y, z;

	if(irank == intExctPortNum){

		// 励振点の設定
		x = intExctLenPart;

		for(y = ex_y_st; y < ex_y_ed; y++){
			for(z = ex_z_st; z < ex_z_ed; z++){
#if _EXITATION_FUNC	// CW励振


				//面内正弦分布励振の場合 ←解析空間が偶数セルか奇数セルかで励振が異なるのでその都度注意

				// スラブ厚の半分のセル数:偶数 導波路幅の半分のセル数:偶数
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt);

				// スラブ厚の半分のセル数:奇数 導波路幅の半分のセル数:奇数
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt);
				Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 01
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st)) * sin(omega0*n*dt); // 02
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); // 03
				//Hz[x][y][z] += cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st - 1)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st - 1)) * sin(omega0*n*dt); // 04
#else	// Gaussian励振

				//Hz[x][(YMAX+1)/2][intSlabCen] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2));
				Hz[x][y][z] += 1000 * cos(omega0*(n-Npeak)*dt) * exp(-(SQ(sigma*dt*(n-Npeak))/2)) * cos(0.5*PI*(y - ex_y_ed + 1)/(ex_y_ed - ex_y_st)) * cos(0.5*PI*(z - ex_z_ed + 1)/(ex_z_ed - ex_z_st));
#endif
			}
		}
	}


	/****************************** 磁界の対称境界条件(4回対称) ******************************/

	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax_ff; z++){
			Hx[x][ymax][z] = Hx[x][ymax-1][z];		// 偶関数
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff+1; z++){
			Hz[x][ymax][z] = Hz[x][ymax-1][z];		// 偶関数
		}
	}
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			Hy[x][y][zmax_ff] = -Hy[x][y][zmax_ff-1];		// 奇関数
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Hx[x][y][zmax_ff] = -Hx[x][y][zmax_ff-1];		// 奇関数
		}
	}

	/****************************** 磁界の対称境界条件(4回対称) ******************************/
}


// モデルへの励振点，観測点の記録
void observation_func(){

	if(irank == intObseInPortNum){ //入射

		for(int x = intObseLenPart1; x < intObseLenPart2; x++){
			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax_ff - intObseHeig; z < zmax_ff; z++){		//矩形導波路断面Z領域判断．
					/****************************** 観測面の修正(2013/8/8) ******************************/
					if((y == YMAX-1) && (z == (intSlabCen-1))){		//活性層断面中央点の時を記録
						//cell[x][y][z] = 4; 					//中央点確認用
					}
					//cell[x][y][z] += OBSERVATION; 		//面確認用
				}
			}
		}
	}
	if (irank == intExctPortNum){
		int x;
		x = intExctLenPart;
		for(int y = ex_y_st; y <= ex_y_ed-1; y++){		//プラス1しているのはセル数の関係
			for(int z = ex_z_st; z <= ex_z_ed-1; z++){
				cell[x][y][z] += EXITATION; 		//励振面確認用
			}
		}
	}

	if(irank == intObseOutPortNum){ //出射 NODE 2
		for(int x = intObseLenPart4; x < intObseLenPart5; x++){
			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax_ff - intObseHeig; z < zmax_ff; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
					//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
					/****************************** 観測面の修正(2013/8/8) ******************************/

					if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
						//cell[x][y][z] = 4; 					// 中央点確認用
					}
					//cell[x][y][z] += OBSERVATION; 			// 面確認用
				}
			}
		}
	}

	if(irank == intObseCenPortNum){ //出射 NODE 2
		int x = intObseLenPart7;
		/****************************** 観測面の修正(2013/8/8) ******************************/
		for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
			for(int z = zmax_ff - intObseHeig; z < zmax_ff; z++){		//矩形導波路断面Z領域判断．
				/****************************** 観測面の修正(2013/8/8) ******************************/

				if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
					//cell[x][y][z] = 4; 					// 中央点確認用
				}
				//cell[x][y][z] += OBSERVATION; 			// 面確認用
			}
		}

	}
}


void calc_efield(){

	double dex, dey, dez;
	double cnstEx, cnstEy, cnstEz;

	// Ex
	for(x = 0; x < xmax; x++){
		for(y = 1; y < ymax+1; y++){		// Exはy軸に対して奇関数
			for(z = 1; z < zmax_ff+1; z++){
				cnstEx = dt / epsilonx[x][y][z];
				dex = ( (Hz[x][y][z] - Hz[x][y-1][z]) / dy) - ( (Hy[x][y][z] - Hy[x][y][z-1]) / dz);
				Ex[x][y][z] = Ex[x][y][z] + cnstEx * dex;
			}
		}
	}

	// Ey
	for(x = 1; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 1; z < zmax_ff+1; z++){
				cnstEy = dt / epsilony[x][y][z];
				dey = ( (Hx[x][y][z] - Hx[x][y][z-1]) / dz)-( (Hz[x][y][z] - Hz[x-1][y][z]) / dx);
				Ey[x][y][z] = Ey[x][y][z] + cnstEy * dey;
			}
		}
	}

	// Ez
	for(x = 1; x < xmax; x++){
		for(y = 1; y < ymax+1; y++){		// Ezはy軸に対して奇関数
			for(z = 0; z < zmax_ff; z++){
				cnstEz = dt / epsilonz[x][y][z];
				dez = ( (Hy[x][y][z] - Hy[x-1][y][z]) / dx) - ( (Hx[x][y][z] - Hx[x][y-1][z]) / dy);
				Ez[x][y][z] = Ez[x][y][z] + cnstEz * dez;
			}
		}
	}


	/****************************** 電界の対称境界条件 ******************************/

	// 境界面で反対称となる電界成分の境界面上の値を0としている
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff+1; z++){
			Ex[x][ymax][z] = 0.0;		// 奇関数
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax_ff; z++){
			Ez[x][ymax][z] = 0.0;		// 奇関数
		}
	}
	/****************************** 電界の対称境界条件 ******************************/
}



void calc_hfield(){

	double dhx, dhy, dhz;

	// Hx
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax_ff; z++){
				dhx = ( (Ey[x][y][z+1] - Ey[x][y][z]) / dz) - ( (Ez[x][y+1][z] - Ez[x][y][z]) / dy);
				Hx[x][y][z] = Hx[x][y][z] + cnstHxyz * dhx;
			}
		}
	}

	// Hy
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax_ff; z++){
				dhy = ( (Ez[x+1][y][z] - Ez[x][y][z]) / dx) - ( (Ex[x][y][z+1] - Ex[x][y][z]) / dz);
				Hy[x][y][z] = Hy[x][y][z] + cnstHxyz * dhy;
			}
		}
	}

	// Hz
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax_ff+1; z++){
				dhz = ((Ex[x][y+1][z] - Ex[x][y][z]) / dy) - ((Ey[x+1][y][z] - Ey[x][y][z]) / dx);
				Hz[x][y][z] = Hz[x][y][z] + cnstHxyz * dhz;
			}
		}
	}
}


// Mur2次，1次の吸収境界条件から端面の計算
void absorpt_bound_condition(){

	/****************************** 対称境界条件 ******************************/

	// 2回対称
	//for(z = 0; z < zmax+1; z++){
	//	Exn1y00[xmax][z] = Exn1y00[xmax-1][z];
	//	Exn1y01[xmax][z] = Exn1y01[xmax-1][z];
	//	Exn1ym0[xmax][z] = Exn1ym0[xmax-1][z];
	//	Exn1ym1[xmax][z] = Exn1ym1[xmax-1][z];
	//}
	//for(y = 0; y < ymax+1; y++){
	//	Exn1z00[xmax][y] = Exn1z00[xmax-1][y];
	//	Exn1z01[xmax][y] = Exn1z01[xmax-1][y];
	//	Exn1zm0[xmax][y] = Exn1zm0[xmax-1][y];
	//	Exn1zm1[xmax][y] = Exn1zm1[xmax-1][y];
	//}
	//for(y = 0; y < ymax; y++){
	//	Eyn1z00[xmax+1][y] = -Eyn1z00[xmax-1][y];
	//	Eyn1z01[xmax+1][y] = -Eyn1z01[xmax-1][y];
	//	Eyn1zm0[xmax+1][y] = -Eyn1zm0[xmax-1][y];
	//	Eyn1zm1[xmax+1][y] = -Eyn1zm1[xmax-1][y];
	//}
	//for(z = 0; z < zmax; z++){
	//	Ezn1y00[xmax+1][z] = -Ezn1y00[xmax-1][z];
	//	Ezn1y01[xmax+1][z] = -Ezn1y01[xmax-1][z];
	//	Ezn1ym0[xmax+1][z] = -Ezn1ym0[xmax-1][z];
	//	Ezn1ym1[xmax+1][z] = -Ezn1ym1[xmax-1][z];
	//}

	// 4回対称
	for(z = 0; z < zmax_ff+1; z++){
		Eyn1x00[ymax][z] = Eyn1x00[ymax-1][z];
		Eyn1x01[ymax][z] = Eyn1x01[ymax-1][z];
		Eyn1xm0[ymax][z] = Eyn1xm0[ymax-1][z];
		Eyn1xm1[ymax][z] = Eyn1xm1[ymax-1][z];
	}
	for(x = 0; x < xmax+1; x++){
		Eyn1z00[x][ymax] = Eyn1z00[x][ymax-1];
		Eyn1z01[x][ymax] = Eyn1z01[x][ymax-1];
		//Eyn1zm0[x][ymax] = Eyn1zm0[x][ymax-1];
		//Eyn1zm1[x][ymax] = Eyn1zm1[x][ymax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1z00[x][ymax+1] = -Exn1z00[x][ymax-1];
		Exn1z01[x][ymax+1] = -Exn1z01[x][ymax-1];
		//Exn1zm0[x][ymax+1] = -Exn1zm0[x][ymax-1];
		//Exn1zm1[x][ymax+1] = -Exn1zm1[x][ymax-1];
	}
	for(z = 0; z < zmax_ff; z++){
		Ezn1x00[ymax+1][z] = -Ezn1x00[ymax-1][z];
		Ezn1x01[ymax+1][z] = -Ezn1x01[ymax-1][z];
		Ezn1xm0[ymax+1][z] = -Ezn1xm0[ymax-1][z];
		Ezn1xm1[ymax+1][z] = -Ezn1xm1[ymax-1][z];
	}

	// 8回対称
	for(y = 0; y < ymax+1; y++){
		Ezn1x00[y][zmax_ff] = -Ezn1x00[y][zmax_ff-1];
		Ezn1x01[y][zmax_ff] = -Ezn1x01[y][zmax_ff-1];
		Ezn1xm0[y][zmax_ff] = -Ezn1xm0[y][zmax_ff-1];
		Ezn1xm1[y][zmax_ff] = -Ezn1xm1[y][zmax_ff-1];
	}
	for(x = 0; x < xmax+1; x++){
		Ezn1y00[x][zmax_ff] = -Ezn1y00[x][zmax_ff-1];
		Ezn1y01[x][zmax_ff] = -Ezn1y01[x][zmax_ff-1];
		//Ezn1ym0[x][zmax] = -Ezn1ym0[x][zmax-1];
		//Ezn1ym1[x][zmax] = -Ezn1ym1[x][zmax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1y00[x][zmax_ff+1] = Exn1y00[x][zmax_ff-1];
		Exn1y01[x][zmax_ff+1] = Exn1y01[x][zmax_ff-1];
		//Exn1ym0[x][zmax+1] = Exn1ym0[x][zmax-1];
		//Exn1ym1[x][zmax+1] = Exn1ym1[x][zmax-1];
	}
	for(y = 0; y <= ymax-1; y++){
		Eyn1x00[y][zmax_ff+1] = Eyn1x00[y][zmax_ff-1];
		Eyn1x01[y][zmax_ff+1] = Eyn1x01[y][zmax_ff-1];
		Eyn1xm0[y][zmax_ff+1] = Eyn1xm0[y][zmax_ff-1];
		Eyn1xm1[y][zmax_ff+1] = Eyn1xm1[y][zmax_ff-1];
	}
	/****************************** 対称境界条件 ******************************/




	/****************************** Murの2次の吸収境界条件(Ex) ******************************/

	double u1ax1, u2ax1,u3ax1, u4ax1;
	double u1bx1, u2bx1,u3bx1, u4bx1;
	double u2xa1;

	double velo_dt;

	if(irank != IRANK_MAX){
		for(x = 1; x < xmax; x++){
			for(z = 1; z < zmax_ff+1; z++){
				velo_dt = (C0 / sqrt(epsilonx[x][0][z]/epsilon0) ) * dt;

				u1ax1 = (velo_dt - dy) / (velo_dt + dy);
				u2ax1 = (2.0 * dy) / (velo_dt + dy);
				u3ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy) );
				u4ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy) );

				Ex[x][0][z] = -Exn2y01[x][z]
				+u1ax1 * (Ex[x][1][z] + Exn2y00[x][z])
					+u2ax1 * (Exn1y00[x][z] + Exn1y01[x][z])
					+u3ax1 * (Exn1y00[x+1][z] - 2.0 * Exn1y00[x][z] + Exn1y00[x-1][z] + Exn1y01[x+1][z] - 2.0 * Exn1y01[x][z] + Exn1y01[x-1][z])
					+u4ax1 * (Exn1y00[x][z+1] - 2.0 * Exn1y00[x][z] + Exn1y00[x][z-1] + Exn1y01[x][z+1] - 2.0 * Exn1y01[x][z] + Exn1y01[x][z-1]);
			}
		}
		for(x = 1; x < xmax; x++){
			for(y = 1; y < ymax+1; y++){
				velo_dt = (C0 / sqrt(epsilonx[x][y][0]/epsilon0) ) * dt;

				u1bx1 = (velo_dt - dz) / (velo_dt + dz);
				u2bx1 = (2.0 * dz) / (velo_dt + dz);
				u3bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) );
				u4bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) );

				Ex[x][y][0] = -Exn2z01[x][y]
				+u1bx1 * (Ex[x][y][1] + Exn2z00[x][y])
					+u2bx1 * (Exn1z00[x][y] + Exn1z01[x][y])
					+u3bx1 * (Exn1z00[x+1][y] - 2.0 * Exn1z00[x][y] + Exn1z00[x-1][y] + Exn1z01[x+1][y] - 2.0*Exn1z01[x][y] + Exn1z01[x-1][y])
					+u4bx1 * (Exn1z00[x][y+1] - 2.0 * Exn1z00[x][y] + Exn1z00[x][y-1] + Exn1z01[x][y+1] - 2.0*Exn1z01[x][y] + Exn1z01[x][y-1]);
			}
		}
	}
	else{
		for(x = 1; x < xmax-1; x++){
			for(z = 1; z < zmax_ff+1; z++){
				velo_dt = (C0 / sqrt(epsilonx[x][0][z]/epsilon0) ) * dt;

				u1ax1 = (velo_dt - dy) / (velo_dt + dy);
				u2ax1 = (2.0 * dy) /(velo_dt + dy);
				u3ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy));
				u4ax1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy));

				Ex[x][0][z] = -Exn2y01[x][z]
				+u1ax1 * (Ex[x][1][z] + Exn2y00[x][z])
					+u2ax1 * (Exn1y00[x][z] + Exn1y01[x][z])
					+u3ax1 * (Exn1y00[x+1][z] - 2.0 * Exn1y00[x][z] + Exn1y00[x-1][z] + Exn1y01[x+1][z] - 2.0 * Exn1y01[x][z] + Exn1y01[x-1][z])
					+u4ax1 * (Exn1y00[x][z+1] - 2.0 * Exn1y00[x][z] + Exn1y00[x][z-1] + Exn1y01[x][z+1] - 2.0 * Exn1y01[x][z] + Exn1y01[x][z-1]);
			}
		}
		for(x = 1; x < xmax-1; x++){
			for(y = 1; y < ymax+1; y++){
				velo_dt = (C0 / sqrt(epsilonx[x][y][0]/epsilon0) ) * dt;

				u1bx1 = (velo_dt - dz) / (velo_dt + dz);
				u2bx1 = (2.0 * dz) / (velo_dt + dz);
				u3bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) );
				u4bx1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) );

				Ex[x][y][0] = -Exn2z01[x][y]
				+u1bx1 * (Ex[x][y][1] + Exn2z00[x][y])
					+u2bx1 * (Exn1z00[x][y] + Exn1z01[x][y])
					+u3bx1 * (Exn1z00[x+1][y] - 2.0 * Exn1z00[x][y] + Exn1z00[x-1][y] + Exn1z01[x+1][y] - 2.0 * Exn1z01[x][y] + Exn1z01[x-1][y])
					+u4bx1 * (Exn1z00[x][y+1] - 2.0 * Exn1z00[x][y] + Exn1z00[x][y-1] + Exn1z01[x][y+1] - 2.0 * Exn1z01[x][y] + Exn1z01[x][y-1]);
			}
		}
	}

	/****************************** Murの2次の吸収境界条件(Ex) ******************************/

	double u1xa, u1xc;
	double u2xa, u2xc;

	/****************************** Murの1次の吸収境界条件(Ex) ******************************/

	for(y = 1; y < ymax+1; y++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonx[0][y][0]/epsilon0) ) * dt;
			u2xa = (velo_dt - dz) / (velo_dt + dz);

			Ex[0][y][0] = Exn1z01[0][y] + u2xa * (Ex[0][y][1] - Exn1z00[0][y]);
		}

		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonx[xmax-1][y][0]/epsilon0) ) * dt;
			u2xc = (velo_dt - dz) / (velo_dt + dz);

			Ex[xmax-1][y][0] = Exn1z01[xmax-1][y] + u2xc * (Ex[xmax-1][y][1] - Exn1z00[xmax-1][y]);
		}
	}

	for(z = 1; z < zmax_ff; z++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonx[0][0][z]/epsilon0) ) * dt;
			u1xa = (velo_dt - dy) / (velo_dt + dy);

			Ex[0][0][z] = Exn1y01[0][z] + u1xa * (Ex[0][1][z] - Exn1y00[0][z]);
		}

		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonx[xmax-1][0][z]/epsilon0) ) * dt;
			u1xc = (velo_dt - dy) / (velo_dt + dy);

			Ex[xmax-1][0][z] = Exn1y01[xmax-1][z] + u1xc * (Ex[xmax-1][1][z] - Exn1y00[xmax-1][z]);
		}
	}


	// 辺(Murの1次の吸収境界条件) -- y平面とz平面からそれぞれ算出される値の平均値を取る
	if (irank != IRANK_MIN){
		for(x = 0; x < xmax; x++){
			velo_dt = (C0 / sqrt(epsilonx[x][0][0]/epsilon0) ) * dt;

			u2xa1 = (velo_dt - dx) / (velo_dt + dx);

			Ex[x][0][0] = 0.5 * (Exn1z01[x][0] + u2xa1 * (Ex[x][0][1] - Exn1z00[x][0])
				+ Exn1y01[x][0] + u2xa1 * (Ex[x][1][0] - Exn1y00[x][0]) );
		}
	}
	else{
		for(x = 1; x < xmax; x++){
			velo_dt = (C0 / sqrt(epsilonx[x][0][0]/epsilon0) ) * dt;

			u2xa1 = (velo_dt - dx) / (velo_dt + dx);

			Ex[x][0][0] = 0.5 * (Exn1z01[x][0] + u2xa1 * (Ex[x][0][1] - Exn1z00[x][0])
				+ Exn1y01[x][0] + u2xa1 * (Ex[x][1][0] - Exn1y00[x][0]) );
		}
	}

	/****************************** Murの1次の吸収境界条件(Ex) ******************************/

	double u1by1, u2by1, u3by1, u4by1;
	double u1cy1, u2cy1, u3cy1, u4cy1;
	double u1cy2, u2cy2, u3cy2, u4cy2;

	/****************************** Murの2次の吸収境界条件(Ey) ******************************/

	for(x = 1; x < xmax; x++){
		for(y = 1; y < ymax; y++){
			velo_dt = (C0 / sqrt(epsilony[x][y][0]/epsilon0) ) * dt;

			u1by1 = (velo_dt - dz) / (velo_dt + dz);
			u2by1 = (2.0 * dz) / (velo_dt + dz);
			u3by1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dz) );
			u4by1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dz) );

			Ey[x][y][0] = -Eyn2z01[x][y]
			+ u1by1 * (Ey[x][y][1] + Eyn2z00[x][y])
				+ u2by1 * (Eyn1z00[x][y] + Eyn1z01[x][y])
				+ u3by1 * (Eyn1z00[x+1][y] - 2.0 * Eyn1z00[x][y] + Eyn1z00[x-1][y] + Eyn1z01[x+1][y] - 2.0 * Eyn1z01[x][y] + Eyn1z01[x-1][y])
				+ u4by1 * (Eyn1z00[x][y+1] - 2.0 * Eyn1z00[x][y] + Eyn1z00[x][y-1] + Eyn1z01[x][y+1] - 2.0 * Eyn1z01[x][y] + Eyn1z01[x][y-1]);
		}
	}
	if(irank == IRANK_MIN){
		for(y = 1; y < ymax; y++){
			for(z = 1; z < zmax_ff+1; z++){
				velo_dt = (C0 / sqrt(epsilony[0][y][z]/epsilon0) ) * dt;

				u1cy1 = (velo_dt - dx) / (velo_dt + dx);
				u2cy1 = (2.0 * dx) / (velo_dt + dx);
				u3cy1 = (dx * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
				u4cy1 = (dx * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

				Ey[0][y][z] = -Eyn2x01[y][z]
				+u1cy1 * (Ey[1][y][z] + Eyn2x00[y][z])
					+u2cy1 * (Eyn1x00[y][z] + Eyn1x01[y][z])
					+u3cy1 * (Eyn1x00[y+1][z] - 2.0 * Eyn1x00[y][z] + Eyn1x00[y-1][z] + Eyn1x01[y+1][z] - 2.0 * Eyn1x01[y][z] + Eyn1x01[y-1][z])
					+u4cy1 * (Eyn1x00[y][z+1] - 2.0 * Eyn1x00[y][z] + Eyn1x00[y][z-1] + Eyn1x01[y][z+1] - 2.0 * Eyn1x01[y][z] + Eyn1x01[y][z-1]);
			}
		}
	}
	if(irank == IRANK_MAX){
		for(y = 1; y < ymax; y++){
			for(z = 1; z < zmax_ff+1; z++){
				velo_dt = (C0 / sqrt(epsilony[xmax][y][z]/epsilon0) ) * dt;

				u1cy2 = (velo_dt - dx) / (velo_dt + dx);
				u2cy2 = (2.0 * dx) / (velo_dt + dx);
				u3cy2 = (dx * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
				u4cy2 = (dx * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

				Ey[xmax][y][z] =  -Eyn2xm1[y][z]
				+ u1cy2 * (Ey[xmax-1][y][z] + Eyn2xm0[y][z])
					+ u2cy2 * (Eyn1xm0[y][z] + Eyn1xm1[y][z])
					+ u3cy2 * (Eyn1xm0[y+1][z] - 2.0 * Eyn1xm0[y][z] + Eyn1xm0[y-1][z] + Eyn1xm1[y+1][z] - 2.0 * Eyn1xm1[y][z] + Eyn1xm1[y-1][z])
					+ u4cy2 * (Eyn1xm0[y][z+1] - 2.0 * Eyn1xm0[y][z] + Eyn1xm0[y][z-1] + Eyn1xm1[y][z+1] - 2.0 * Eyn1xm1[y][z] + Eyn1xm1[y][z-1]);
			}
		}
	}

	/****************************** Murの2次の吸収境界条件(Ey) ******************************/


	double u2ya, u3ya, u3yb;
	double u2ya1;
	double u2yc1;

	/****************************** Murの1次の吸収境界条件(Ey) ******************************/

	for(x = 1; x < xmax; x++){
		velo_dt = (C0 / sqrt(epsilony[x][0][0]/epsilon0) ) * dt;
		u2ya = (velo_dt - dz) / (velo_dt + dz);
		Ey[x][0][0] = Eyn1z01[x][0] + u2ya * (Ey[x][0][1] - Eyn1z00[x][0]);
	}

	for(z = 1; z < zmax_ff+1; z++){
		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilony[0][0][z]/epsilon0) ) * dt;

			u3ya = (velo_dt - dx) / (velo_dt + dx);
			Ey[0][0][z] = Eyn1x01[0][z] + u3ya * (Ey[1][0][z] - Eyn1x00[0][z]);
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilony[xmax][0][0]/epsilon0) ) * dt;

			u3yb = (velo_dt - dx) / (velo_dt + dx);
			Ey[xmax][0][z] = Eyn1xm1[0][z] + u3yb * (Ey[xmax-1][0][z] - Eyn1xm0[0][z]);
		}
	}

	// 辺(Murの1次の吸収境界条件) --x平面とz平面からそれぞれ算出される値の平均値を取る
	for(y = 0; y < ymax; y++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonx[0][y][0]/epsilon0) ) * dt;

			u2ya1 = (velo_dt - dz) / (velo_dt + dz);
			Ey[0][y][0] = 0.5 * (Eyn1z01[0][y] + u2ya1 * (Ey[0][y][1] - Eyn1z00[0][y])
				+ Eyn1x01[y][0] + u2ya1 * (Ey[1][y][0] - Eyn1x00[y][0]));
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilony[xmax][y][0]/epsilon0) ) * dt;

			u2yc1 = (velo_dt - dz) / (velo_dt + dz);
			Ey[xmax][y][0] = 0.5*(Eyn1z01[xmax][y] + u2yc1 * (Ey[xmax][y][1] - Eyn1z00[xmax][y])
				+ Eyn1xm1[y][0] + u2yc1 * (Ey[xmax-1][y][0] - Eyn1xm0[y][0]));
		}
	}
	/****************************** Murの1次の吸収境界条件(Ey) ******************************/

	double u1az1, u2az1, u3az1, u4az1;
	double u1cz1, u2cz1, u3cz1, u4cz1;
	double u1cz2, u2cz2, u3cz2, u4cz2;

	/****************************** Murの2次の吸収境界条件(Ez) ******************************/

	for(x = 1; x < xmax; x++){
		for(z = 1; z < zmax_ff; z++){
			velo_dt = (C0 / sqrt(epsilonz[x][0][z] / epsilon0) ) * dt;

			u1az1 = (velo_dt - dy) / (velo_dt + dy);
			u2az1 = (2.0 * dy) / (velo_dt + dy);
			u3az1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dx) * (velo_dt + dy) );
			u4az1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dy) );

			Ez[x][0][z] = -Ezn2y01[x][z]
			+ u1az1 * (Ez[x][1][z] + Ezn2y00[x][z])
				+ u2az1 * (Ezn1y00[x][z] + Ezn1y01[x][z])
				+ u3az1 * (Ezn1y00[x+1][z] - 2.0 * Ezn1y00[x][z] + Ezn1y00[x-1][z] + Ezn1y01[x+1][z] - 2.0 * Ezn1y01[x][z] + Ezn1y01[x-1][z])
				+ u4az1 * (Ezn1y00[x][z+1] - 2.0 * Ezn1y00[x][z] + Ezn1y00[x][z-1] + Ezn1y01[x][z+1] - 2.0 * Ezn1y01[x][z] + Ezn1y01[x][z-1]);
		}
	}

	for(y = 1; y < ymax+1; y++){
		for(z = 1; z < zmax_ff; z++){
			if(irank == IRANK_MIN){
				velo_dt = (C0 / sqrt(epsilonz[0][y][z] / epsilon0) ) * dt;

				u1cz1 = (velo_dt - dx) / (velo_dt + dx);
				u2cz1 = (2.0 * dx) / (velo_dt + dx);
				u3cz1 = (dy * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
				u4cz1 = (dz * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

				Ez[0][y][z] = -Ezn2x01[y][z]
				+ u1cz1 * (Ezn2x00[y][z] + Ez[1][y][z])
					+ u2cz1 * (Ezn1x00[y][z] + Ezn1x01[y][z])
					+ u3cz1 * (Ezn1x00[y+1][z] - 2.0 * Ezn1x00[y][z] + Ezn1x00[y-1][z] + Ezn1x01[y+1][z] - 2.0 * Ezn1x01[y][z] + Ezn1x01[y-1][z])
					+ u4cz1 * (Ezn1x00[y][z+1] - 2.0 * Ezn1x00[y][z] + Ezn1x00[y][z-1] + Ezn1x01[y][z+1] - 2.0 * Ezn1x01[y][z] + Ezn1x01[y][z-1]);
			}
			if(irank == IRANK_MAX){
				velo_dt = (C0 / sqrt(epsilonz[xmax][y][z] / epsilon0) ) * dt;

				u1cz2 = (velo_dt - dx) / (velo_dt + dx);
				u2cz2 = (2.0 * dx) / (velo_dt + dx);
				u3cz2 = (dy * SQ (velo_dt)) / (2.0 * SQ(dy) * (velo_dt + dx) );
				u4cz2 = (dz * SQ (velo_dt)) / (2.0 * SQ(dz) * (velo_dt + dx) );

				Ez[xmax][y][z] = -Ezn2xm1[y][z]
				+ u1cz2 * (Ezn2xm0[y][z] + Ez[xmax-1][y][z])
					+ u2cz2 * (Ezn1xm1[y][z] + Ezn1xm0[y][z])
					+ u3cz2 * (Ezn1xm1[y+1][z] - 2.0 * Ezn1xm1[y][z] + Ezn1xm1[y-1][z] + Ezn1xm0[y+1][z] - 2.0 * Ezn1xm0[y][z] + Ezn1xm0[y-1][z])
					+ u4cz2 * (Ezn1xm1[y][z+1] - 2.0 * Ezn1xm1[y][z] + Ezn1xm1[y][z-1] + Ezn1xm0[y][z+1] - 2.0 * Ezn1xm0[y][z] + Ezn1xm0[y][z-1]);
			}
		}
	}

	/****************************** Murの2次の吸収境界条件(Ez) ******************************/

	double u1za, u3za, u3zb;
	double u1za1, u1zb1;

	/****************************** Murの1次の吸収境界条件(Ez) ******************************/
	for(x = 1; x < xmax; x++){
		velo_dt = (C0 / sqrt(epsilonz[x][0][0] / epsilon0) ) * dt;
		u1za = (velo_dt - dy) / (velo_dt + dy);

		Ez[x][0][0] = Ezn1y01[x][0] + u1za * (Ez[x][1][0] - Ezn1y00[x][0]);
	}

	for(y = 1; y < ymax+1; y++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonz[0][y][0] / epsilon0) ) * dt;
			u3za = (velo_dt - dx) / (velo_dt + dx);

			Ez[0][y][0] = Ezn1x01[y][0] + u3za * (Ez[1][y][0] - Ezn1x00[y][0]);
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonz[xmax][y][0] / epsilon0) ) * dt;
			u3zb = (velo_dt - dx) / (velo_dt + dx);

			Ez[xmax][y][0] = Ezn1xm1[y][0] + u3zb * (Ez[xmax-1][y][0] - Ezn1xm0[y][0]);
		}
	}

	// 辺(Murの1次の吸収境界条件) --x平面とy平面からそれぞれ算出される値の平均値を取る
	for(z = 0; z < zmax_ff+1; z++){

		if(irank == IRANK_MIN){
			velo_dt = (C0 / sqrt(epsilonz[0][0][z] / epsilon0) ) * dt;
			u1za1 = (velo_dt - dy) / (velo_dt + dy);

			Ez[0][0][z] = 0.5 * (Ezn1y01[0][z] + u1za1 * (Ez[0][1][z] - Ezn1y00[0][z])
				+ Ezn1x01[0][z] + u1za1 * (Ez[1][0][z] - Ezn1x00[0][z]) );
		}
		if(irank == IRANK_MAX){
			velo_dt = (C0 / sqrt(epsilonz[xmax][0][z] / epsilon0) ) * dt;
			u1zb1 = (velo_dt - dy) / (velo_dt + dy);

			Ez[xmax][0][z] = 0.5 * (Ezn1y01[xmax][z] + u1zb1 * (Ez[xmax][1][z] - Ezn1y00[xmax][z])
				+ Ezn1xm1[0][z] + u1zb1 * (Ez[xmax-1][0][z] - Ezn1xm0[0][z]));
		}
	}
	/****************************** Murの1次の吸収境界条件(Ez) ******************************/

}



/*電界の保存*/
void saving_electric_field(){

	// Ex
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax_ff+1; z++){
			Exn2y00[x][z] = Exn1y00[x][z];
			Exn1y00[x][z] = Ex[x][0][z];
			Exn2y01[x][z] = Exn1y01[x][z];
			Exn1y01[x][z] = Ex[x][1][z];
			//exn2ym1[x][z] = exn1ym1[x][z];
			//exn1ym1[x][z] = Ex[x][ymax-1][z];
			//exn2ym0[x][z] = exn1ym0[x][z];
			//exn1ym0[x][z] = Ex[x][ymax][z];
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			Exn2z00[x][y] = Exn1z00[x][y];
			Exn1z00[x][y] = Ex[x][y][0];
			Exn2z01[x][y] = Exn1z01[x][y];
			Exn1z01[x][y] = Ex[x][y][1];
		}
	}

	// Ey
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Eyn2z00[x][y] = Eyn1z00[x][y];
			Eyn1z00[x][y] = Ey[x][y][0];
			Eyn2z01[x][y] = Eyn1z01[x][y];
			Eyn1z01[x][y] = Ey[x][y][1];
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax_ff+1; z++){
			if(irank == IRANK_MIN){
				Eyn2x00[y][z] = Eyn1x00[y][z];
				Eyn1x00[y][z] = Ey[0][y][z];
				Eyn2x01[y][z] = Eyn1x01[y][z];
				Eyn1x01[y][z] = Ey[1][y][z];
			}
			if(irank == IRANK_MAX){
				Eyn2xm1[y][z] = Eyn1xm1[y][z];
				Eyn1xm1[y][z] = Ey[xmax-1][y][z];
				Eyn2xm0[y][z] = Eyn1xm0[y][z];
				Eyn1xm0[y][z] = Ey[xmax][y][z];
			}
		}
	}

	//Ez
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax_ff; z++){
			Ezn2y00[x][z] = Ezn1y00[x][z];
			Ezn1y00[x][z] = Ez[x][0][z];
			Ezn2y01[x][z] = Ezn1y01[x][z];
			Ezn1y01[x][z] = Ez[x][1][z];
		}
	}
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax_ff; z++){
			if(irank == IRANK_MIN){
				Ezn2x00[y][z] = Ezn1x00[y][z];
				Ezn1x00[y][z] = Ez[0][y][z];
				Ezn2x01[y][z] = Ezn1x01[y][z];
				Ezn1x01[y][z] = Ez[1][y][z];
			}
			if(irank == IRANK_MAX){
				Ezn2xm1[y][z] = Ezn1xm1[y][z];
				Ezn1xm1[y][z] = Ez[xmax-1][y][z];
				Ezn2xm0[y][z] = Ezn1xm0[y][z];
				Ezn1xm0[y][z] = Ez[xmax][y][z];
			}
		}
	}
}


//モデルの出力
void output_model(){

	int tag2 = 2;
	int x, y, z;

#if _FDTD
	/****************************** 計算実行時 ******************************/
	int node;

	MPI_Status status;

	z = intSlabCen - 1;

	// XY平面
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			cell_xy[x][y] = cell[x][y][z];
		}
	}

	if(irank == IRANK_MIN){
		for(x = 0; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(allmodel_xy, "%d\t", cell_xy[x][y]);
				fprintf(model_xy, "%d\t", cell_xy[x][y]);
			}
			fprintf(allmodel_xy, "\n");
			fprintf(model_xy, "\n");
		}
		fclose(model_xy);
	}

	// XZ平面
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff; z++){
			cell_xz[x][z] = cell[x][y][z];
		}
	}

	if(irank == IRANK_MIN){
		for(x = 0; x < xmax; x++){
			for(z = 0; z < zmax_ff; z++){
				fprintf(allmodel_xz, "%d\t", cell_xz[x][z]);
				fprintf(model_xz, "%d\t", cell_xz[x][z]);
			}
			fprintf(allmodel_xz, "\n");
			fprintf(model_xz, "\n");
		}
		fclose(model_xz);
	}

	// YZ平面
	if(irank == NODE/2){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax_ff; z++){
				cell_yz[y][z] = cell[10][y][z];
				fprintf(model_yz, "%d\t", cell_yz[y][z]);
			}
			fprintf(model_yz, "\n");
		}
		fclose(model_yz);
	}


	// それぞれ分割部のモデル
	if(irank != IRANK_MIN){
		MPI_Send(&cell_xy[0][0], (xmax)*(ymax), MPI_INT, 0, tag2, MPI_COMM_WORLD);
		for(x = 1; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(model_xy, "%d\t", cell_xy[x][y]);
			}
			fprintf(model_xy, "\n");
		}
		fclose(model_xy);
	}


	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){
			if(node == IRANK_MAX){

				// 最終段は"のりしろ"が無いので，x方向を-1する
				MPI_Recv(&cell_xy[0][0], (xmax-1)*(ymax), MPI_INT, node, tag2, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax-1; x++){
					for(y = 0; y < ymax; y++){
						fprintf(allmodel_xy, "%d\t", cell_xy[x][y]);
					}
					fprintf(allmodel_xy, "\n");
				}
			}
			else{

				MPI_Recv(&cell_xy[0][0], (xmax)*(ymax), MPI_INT, node, tag2, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax; y++){
						fprintf(allmodel_xy, "%d\t", cell_xy[x][y]);
					}
					fprintf(allmodel_xy, "\n");
				}
			}
		}
		fclose(allmodel_xy);
	}




	if(irank != IRANK_MIN){
		MPI_Send(&cell_xz[0][0], (xmax)*(zmax_ff), MPI_INT, 0, tag2, MPI_COMM_WORLD);
		for(x = 1; x < xmax; x++){
			for(z = 0; z < zmax_ff; z++){
				fprintf(model_xz, "%d\t", cell_xz[x][z]);
			}
			fprintf(model_xz, "\n");
		}
		fclose(model_xz);
	}

	// 全体モデル生成
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){
			if(node == IRANK_MAX){
				// 最終段は"のりしろ"が無いので，x方向を-1する
				MPI_Recv(&cell_xz[0][0], (xmax-1)*(zmax_ff), MPI_INT, node, tag2, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax-1; x++){
					for(z = 0; z < zmax_ff; z++){
						fprintf(allmodel_xz, "%d\t", cell_xz[x][z]);
					}
					fprintf(allmodel_xz, "\n");
				}
			}
			else{
				MPI_Recv(&cell_xz[0][0], (xmax)*(zmax_ff), MPI_INT, node, tag2, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax; x++){
					for(z = 0; z < zmax_ff; z++){
						fprintf(allmodel_xz, "%d\t", cell_xz[x][z]);
					}
					fprintf(allmodel_xz, "\n");
				}
			}
		}
		fclose(allmodel_xz);
	}



	/****************************** 計算実行時 ******************************/

#else
	/****************************** モデル確認時 ******************************/
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			cell_xy[x][y] = cell[x][y][intSlabCen-1];
		}
	}

	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			fprintf(allmodel_xy, "%d\t", cell_xy[x][y]);
			fprintf(model_xy, "%d\t", cell_xy[x][y]);
		}
		fprintf(allmodel_xy, "\n");
		fprintf(model_xy, "\n");
	}
	fclose(model_xy);

	x = intObseLenPart1;
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax; z++){
			cell_yz[y][z] = cell[x][y][z];
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax; z++){
			fprintf(allmodel_yz1, "%d\t", cell_yz[y][z]);
		}
		fprintf(allmodel_yz1, "\n");
	}
	fclose(allmodel_yz1);

	if(irank == intObseOutPortNum){ // 出射
		x = intObseLenPart4;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z];
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz4, "%d\t", cell_yz[y][z]);
			}
			fprintf(allmodel_yz4, "\n");
		}
		fclose(allmodel_yz4);
	}

	if(irank == intObseCenPortNum){			// 中央
		x = intObseLenPart7;
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				cell_yz[y][z] = cell[x][y][z];
			}
		}
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax; z++){
				fprintf(allmodel_yz7, "%d\t", cell_yz[y][z]);
			}
			fprintf(allmodel_yz7, "\n");
		}
		fclose(allmodel_yz7);
	}
	/****************************** モデル確認時 ******************************/
#endif

}


void output_field_write(char *dir_name_def){

	char fname1[40], fname2[40], fname3[40], dir_name[50]; 	//ファイル名格納変数
	int node;
	int tag3 = 3;
	int pi1, pj1, pk1;
	MPI_Status status;
	FILE *HZ1, *HZ1_NODE;
	FILE *HZ2, *HZ2_NODE;
	FILE *HZ3;
	//FILE *HZ1, *HZ2;

	pi1 = x_cen;
	pj1 = y_cen;
	pk1 = z_cen;

	printf("n = %d\n", n);

	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			field_xy[x][y] = Hz[x][y][ex_z_ed-1]; 	//全てのノードで電磁界成分を2次元配列に格納する．
		}
	}


	// モデル出力ファイルポインタの初期化
	if(irank == IRANK_MIN){
		sprintf(fname1, "/Field_Hz_XY_%d_01.txt", n);
		HZ1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname1), "w");
		for(x = 0; x < xmax; x++){
			for(y = 0; y < ymax; y++){
				fprintf(HZ1, "%e\t", field_xy[x][y]);
			}
			fprintf(HZ1, "\n");
		}
	}
	// モデルをホストに送信
	else{
		if(irank != IRANK_MAX){
			MPI_Send(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 		// ノード0以外のノードがノード0に電磁界成分を送る．
		}
		if(irank == IRANK_MAX){
			MPI_Send(&field_xy[0][0], (xmax-1)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// ノード0以外のノードがノード0に電磁界成分を送る．
		}
	}

	// 各ノードごとにXY平面の磁界分布の出力
/*	sprintf(fname1, "/Field_Hz_XY_%d_%d_01.txt", irank, n);
	HZ1_NODE = fopen(strcat(strcpy(dir_name, dir_name_def), fname1), "w");
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			fprintf(HZ1_NODE, "%e\t", field_xy[x][y]);
		}
		fprintf(HZ1_NODE, "\n");
	}*/



	// 受信したモデルから全モデルを作成
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){		// ノード0がノード1から順にデータを受け取り出力していく．
			if(node == IRANK_MAX){					// ノードisize-1のみ1セル小さく設定しているため条件文で分岐
				MPI_Recv(&field_xy[0][0], (xmax-1)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax-1; x++){
					for(y = 0; y < ymax; y++){
						fprintf(HZ1, "%e\t", field_xy[x][y]);
					}
					fprintf(HZ1, "\n");
				}
			}
			else{
				MPI_Recv(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax; x++){
					for(y = 0; y < ymax; y++){
						fprintf(HZ1, "%e\t", field_xy[x][y]);
					}
					fprintf(HZ1, "\n");
				}

			}
		}
		// ファイルポインタを閉じる
		fclose(HZ1);
	}


	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff; z++){
			field_xz[x][z] = Hz[x][YMAX][z]; 	//全てのノードで電磁界成分を2次元配列に格納する．
		}
	}

	// モデル出力ファイルポインタの初期化
	if(irank == IRANK_MIN){
		sprintf(fname2, "/Field_Hz_XZ_%d_01.txt", n);
		HZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname2), "w");
		for(x = 0; x < xmax; x++){
			for(z = 0; z < zmax_ff; z++){
				fprintf(HZ2, "%e\t", field_xz[x][z]);
			}
			fprintf(HZ2, "\n");
		}
	}else{
		if(irank != IRANK_MAX){
			MPI_Send(&field_xz[0][0], (xmax)*(zmax_ff), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD);		// ノード0以外のノードがノード0に電磁界成分を送る．
		}
		if(irank == IRANK_MAX){
			MPI_Send(&field_xz[0][0], (xmax-1)*(zmax_ff), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// ノード0以外のノードがノード0に電磁界成分を送る．
		}
	}

	// 各ノードごとにXZ平面の磁界分布の出力
/*	sprintf(fname2, "/Field_Hz_XZ_%d_%d_01.txt", irank, n);
	HZ2_NODE = fopen(strcat(strcpy(dir_name, dir_name_def), fname2), "w");
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax_ff; z++){
			fprintf(HZ2_NODE, "%e\t", field_xz[x][z]);
		}
		fprintf(HZ2_NODE, "\n");
	}*/

	// 受信したモデルから全モデルを作成
	if(irank == IRANK_MIN){
		for(node = 1; node < ISIZE; node++){		// ノード0がノード1から順にデータを受け取り出力していく．
			if(node == IRANK_MAX){					// ノードisize-1のみ1セル小さく設定しているため条件文で分岐
				MPI_Recv(&field_xz[0][0], (xmax-1)*(zmax_ff), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax-1; x++){
					for(z = 0; z < zmax_ff; z++){
						fprintf(HZ2, "%e\t", field_xz[x][z]);
					}
					fprintf(HZ2, "\n");
				}
			}
			else{
				MPI_Recv(&field_xz[0][0], (xmax)*(zmax_ff), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for(x = 1; x < xmax; x++){
					for(z = 0; z < zmax_ff; z++){
						fprintf(HZ2, "%e\t", field_xz[x][z]);
					}
					fprintf(HZ2, "\n");
				}
			}
		}
		// ファイルポインタを閉じる
		fclose(HZ2);
	}


	if(irank == NODE/2){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax_ff; z++){
				field_yz[y][z] = Hz[10][y][z]; 	//全てのノードで電磁界成分を2次元配列に格納する．
			}
		}

		sprintf(fname3, "/Field_Hz_YZ_%d_01.txt", n);
		HZ3 = fopen(strcat(strcpy(dir_name, dir_name_def), fname3), "w");
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax_ff; z++){
				fprintf(HZ3, "%e\t", field_yz[y][z]);
			}
			fprintf(HZ3, "\n");
		}
	}


	// YZ平面の電界分布の出力
	/*int x;
	double E_yz;
	FILE *EYZ1, *EYZ2, *EYZ3;
	char fname2[40], fname3[40], fname4[40];*/

	/*if(irank == intObseInPortNum){ //入射
		x = intObseLenPart1;
		sprintf(fname2, "/Field_E_YZ_%d_01.txt", n);
		EYZ1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname2), "w");

		for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
			for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ1, "%e\t", E_yz);
			}
			fprintf(EYZ1, "\n");
		}
		fclose(EYZ1);
	}

	if(irank == intObseOutPortNum){			// 出射
		x = intObseLenPart4;
		sprintf(fname3, "/Field_E_YZ_%d_04.txt", n);
		EYZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname3), "w");

		for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
			for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ2, "%e\t", E_yz);
			}
			fprintf(EYZ2, "\n");
		}
		fclose(EYZ2);
	}

	if(irank == intObseCenPortNum){			// 出射
		x = intObseLenPart7;
		sprintf(fname4, "/Field_E_YZ_%d_07.txt", n);
		EYZ3 = fopen(strcat(strcpy(dir_name, dir_name_def), fname4), "w");

		for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
			for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ3, "%e\t", E_yz);
			}
			fprintf(EYZ3, "\n");
		}
		fclose(EYZ3);
	}*/

}

//ファイル出力
void output_field(char *dir_name_def){

	if(n <= Nmax - Fcut){
		// 動作確認のためのファイル出力
		if(n == Ncheck){
			output_field_write (dir_name_def);
		}

		// 定期的なファイル出力
		if(n % Ncutfield == 0){
			output_field_write (dir_name_def);
		}
	}
	if((n >= Nmax - Fcut) && (n <= Nmax)){

		// 安定点でのファイル出力
		if(n % Ncutfield2 == 0){
			output_field_write (dir_name_def);
		}
	}
}


void mcircle(int x_circ, int y_circ, int z_circ, int type){

	double R;
	double Rs,Rb;

	Rs = ((dblRadius_s*1.0e10)/(dx*1.0e10));
	Rb = ((dblRadius_b*1.0e10)/(dx*1.0e10));

	//半径セル数の計算
	if(type == 1)	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています
	else if(type == 2)	R = ((dblRadius2*1.0e10)/(dx*1.0e10));
	else if(type == 3)	R = ((dblRadius3*1.0e10)/(dx*1.0e10));
	else if(type == 5)	R = ((dblRadius5*1.0e10)/(dx*1.0e10));
	else if(type == 6)	R = ((dblRadius6*1.0e10)/(dx*1.0e10));
	else if(type == 7)	R = ((dblRadius7*1.0e10)/(dx*1.0e10));
	else if(type == 8)	R = ((dblRadius8*1.0e10)/(dx*1.0e10));
	else			R = ((dblRadius4*1.0e10)/(dx*1.0e10));



	if(flag_2r == 0){
		rightquartercircle1(x_circ, y_circ, z_circ, type, R);
		leftquartercircle1(x_circ-1, y_circ, z_circ, type, R);
		rightquartercircle2(x_circ, y_circ-1, z_circ, type, R);
		leftquartercircle2(x_circ-1, y_circ-1, z_circ, type, R);
	}else if(flag_2r == 2){
		rightquartercircle1(x_circ, y_circ, z_circ, type, Rs);
		leftquartercircle1(x_circ-1, y_circ, z_circ, type, Rs);
		rightquartercircle2(x_circ, y_circ-1, z_circ, type, Rs);
		leftquartercircle2(x_circ-1, y_circ-1, z_circ, type, Rs);
	}else if(flag_2r == 1){
		rightquartercircle1(x_circ, y_circ, z_circ, type, Rb);
		leftquartercircle1(x_circ-1, y_circ, z_circ, type, Rb);
		rightquartercircle2(x_circ, y_circ-1, z_circ, type, Rb);
		leftquartercircle2(x_circ-1, y_circ-1, z_circ, type, Rb);
	}
}


void halfcircle(int x_circ, int y_circ, int z_circ, int type){

	double R;

	//半径セル数の計算
	if(type == 1)	R = ((dblRadius*1.0e10)/(dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています
	else if(type == 2)	R = ((dblRadius2*1.0e10)/(dx*1.0e10));
	else if(type == 3)	R = ((dblRadius3*1.0e10)/(dx*1.0e10));
	else if(type == 5)	R = ((dblRadius5*1.0e10)/(dx*1.0e10));
	else if(type == 6)	R = ((dblRadius6*1.0e10)/(dx*1.0e10));
	else if(type == 7)	R = ((dblRadius7*1.0e10)/(dx*1.0e10));
	else if(type == 8)	R = ((dblRadius8*1.0e10)/(dx*1.0e10));
	else			R = ((dblRadius4*1.0e10)/(dx*1.0e10));

	rightquartercircle2(x_circ, y_circ-1, z_circ, type, R);
	leftquartercircle2(x_circ-1, y_circ-1, z_circ, type, R);
}


void rightquartercircle1(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ+R-1);
	Je = (int) (y_circ+R-1);
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y <= Je; y++){
			r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ+1) * (y-y_circ+1)) ) - 0.5;
			if(r <= R){
				if(type == 1 || type == 3|| type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}


void leftquartercircle1(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ-R+1);
	Je = (int) (y_circ+R-1);
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y <= Je; y++){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ+1) * (y-y_circ+1)) ) - 0.5;
			if(r <= R){
				if(type == 1 || type == 3|| type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}


void rightquartercircle2(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ+R-1);
	Je = (int) (y_circ-R+1);
	for(x = x_circ; x <= Ie; x++){
		for(y = y_circ; y >= Je; y--){
			r = sqrt(double((x-x_circ+1) * (x-x_circ+1) + (y-y_circ-1) * (y-y_circ-1)) ) - 0.5;
			if(r <= R){
				if(type == 1 || type == 3|| type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}

void leftquartercircle2(int x_circ, int y_circ, int z_circ, int type, double R){

	int x, y, Ie, Je;
	double r;

	Ie = (int) (x_circ-R+1);
	Je = (int) (y_circ-R+1);
	for(x = x_circ; x >= Ie; x--){
		for(y = y_circ; y >= Je; y--){
			r = sqrt(double((x-x_circ-1) * (x-x_circ-1) + (y-y_circ-1) * (y-y_circ-1)) ) - 0.5;
			if(r <= R){
				if(type == 1 || type == 3 || type == 5 || type == 6 || type == 7 || type == 8){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else if(type == 2){
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}else{
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}
