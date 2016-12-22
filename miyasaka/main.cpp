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

//PCW_WIDTH_CHIRP
//PITCH_SHIFT_PER
//intPitchShiftPcwPer

#define _FDTD 1// FDTD計算			0 : モデル掃き出し(プリプロセッサでコンパイルを変更させる)
//									1 : 計算実行

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

#define BOUNDARYLINE 31//SiWGとPWWの境界条件 0 : 従来WG構造 2015/11/11☆
//											 1 : WGテーパ円孔1列目(線形) 2015/11/19
//	　　　　							     2 : WGテーパ円孔1列目(放物線)
//										     3 : WGテーパ円孔1列目(逆放物線)2015/11/24
//										     4 : WGテーパ円孔1列目(放物線改)2015/11/25　 出射のテーパ緩やか　未計算
//										     5 : WGテーパ円孔1列目(線形改)2015/11/25　　
//										     6 : WGテーパ円孔1列目(放物線改2)2015/11/25　入出射のテーパ緩やかに
//										     7 : WGテーパ円孔1列目(放物線改3)2015/11/26　入出射を小さく
//										     8 : テーパ2015/12/7 國分先生 12/14これは保留
//										     9 : WG緩やかテーパ(出射のみ)2015/12/14　
//										    10 : さらにWG緩やかテーパ(出射のみ)2015/12/14
//										    11 : PCW近接オフセット兼さらにWG緩やかテーパ(出射のみ)2015/12/15
//										    12 : PCW近接オフセット深く兼さらにWG緩やかテーパ(出射のみ)2015/12/15
//										    13 : 12+入射放物線2015/12/15
//										    14 : 12+モニター後ろ
//										    141 : モニター後ろのみ 2015/12/16
//										    15 : 14+テーパ変化 2015/12/16
//										    16 : 14+テーパ変化2 2015/12/16
//										    17 : 14+テーパ変化3 2015/12/24 傾き0.2
//										    18 : 14+テーパ変化4 2015/12/24 傾き0.12
//										    19 : 14+テーパ変化5 2015/12/25 傾き0.10
//										    20 : 14+テーパ変化6 2015/12/25 傾き0.08
//										    21 : 20 + 出射放物線テーパ 2016/1/4 傾き0.08
//										    22 : 14+テーパ変化7 2016/1/6 傾き0.12
//										    23 : 14+テーパ変化8 2016/1/6 傾き0.16
//										    24 : 14+テーパ変化9 2016/1/6 傾き0.20
//										    25 : 14+テーパ変化10 2016/1/7 傾き0.24
//										    26 : 14+テーパ変化11 2016/1/7 傾き0.04
//										    27 : 14+テーパ変化12 2016/1/7 傾き0.02 実際はノーテーパー(テーパは生じないので)
//										    28 : ノーテーパー　2016/2/5y方向長尺(Ymax=155，SY =-231にする必要あり) 従来より10セル分pcw端のスペース(14)
//										    29 : ノーテーパー　2016/2/5y方向長尺(Ymax=157，SY =-273にする必要あり) 従来より12セル分pcw端のスペース(16)
//										    30 : 2016/2/8 (これは容量が大きいのでモデルを見るにはXmaxを減らす必要あり)ノーテーパー　y方向長尺(Ymax=162，SY =-378にする必要あり) 従来より17セル分pcw端のスペース(21)
//										    31 : 2016/2/8 (?)ノーテーパー y方向長尺(Ymax=172，SY =-588)にする必要あり) 従来より27セル分pcw端のスペース(31)　☆これだと横の分布が良く見える
//										    32 : 2016/2/16 31+ 出射微小テーパ
//										    33 : 2016/3/3 31+ 入射微小テーパ(長尺導波路)
//34 : 2016/10/11入出射円錐状構造 dtaper = 168 nm, Ymax = 172，SY = -588に固定 傾き0.08 WIRE_WID_OFFSET, OUT = 315 (24→31 nmに約30%広げる) N = 7 テーパ0.08
//35 : 34+ テーパ0.06
//36 : 34+ テーパ0.04
//									    37 : 2016/10/18 斜めPCW入射
//38 : 2016/10/18 斜めPCW入射 幅テーパ(168 nm)
//39 : 2016/10/20 斜めPCW入出射
//40 : 2016/10/20 斜めPCW入出射 より低傾き
//41 : 2016/10/21 斜めPCW入出(39) 円孔位置調節 外側接触(y方向2.5周期) 論文と同様の構造
//42 : 2016/10/21 斜めPCW入出(39) 円孔位置調節 中心接触(y方向2.5周期)
//43 : 2016/10/21 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2.5周期)
//44 : 2016/10/24 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2周期)
//45 : 2016/10/25 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向1.5周期)
//46 : 2016/10/25 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向1周期)未作成→保留
//47 : 2016/10/25 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2.5周期) PCW周期数　x方向8周期 (41～46は7周期)
//48 : 2016/10/25 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2.5周期) PCW周期数　x方向9周期 (41～46は7周期)
//49 : 2016/11/4 細線評価用  ノーテーパー(31) + PCW_WID 1
//50 : 2016/11/9 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2.5周期) 幅チャ(168 nm)併用 N = 7
//51 : 2016/11/9 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2.5周期) 幅チャ(168 nm)併用 N = 9
//52 : 2016/11/9 斜めPCW入出(41) 円孔位置調節 外側接触 引き離し(y方向2.5周期) 幅チャ(168 nm)併用 N = 9 幅チャ途中から

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
//void magnitude();
void output_field(char*);
void output_field_write(char *);
void output_model();
void calc_energy();
void calc_power();
void calc_powerHz();
void calc_poynting_power();
void calc_poynting_powerHz();

int main(int argc, char **argv){

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION
	double s_time, e_time;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	char time[9];
	char stime[9];
	char etime[9];
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
			//fprintf("Start Time:\t %s\n", time);

			//stime[0] == 0;//☆計算開始時刻をホールド　動作しない
			/*stime[1] == 0;
			stime[2] == 0;
			stime[3] == 0;
			stime[4] == 0;
			stime[5] == 0;
			stime[6] == 0;
			stime[7] == 0;
			stime[8] == 0;
			stime[9] == 0;*/ //1行空けるはず

			//s_time = MPI_Wtime();
			//fprintf(fpparameter, "Start Time1:\t %s\n", s_time);
		}
		// 電磁界計算
		for(n = 1 ; n <= Nmax; n++){

			// 時間ステップ数の表示
			if(n % Ncut == 0){
				//_strtime(time);
				printf("n = %d, \t\t", n);
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

			MPI_Sendrecv( &Ex[1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_send,
				&Ex[xmax][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ey[1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_send,
				&Ey[xmax][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ez[1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_send,
				&Ez[xmax][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ex[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
				&Ex[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ey[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
				&Ey[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Ez[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
				&Ez[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

			int var_i, var_j, var_k;
			for(var_i = 0; var_i <= xmax; var_i++) {
				for(var_j = 0; var_j <= ymax; var_j++) {
					for(var_k = 0; var_k <= zmax; var_k++) {
						printf("(%e %e %e)", Ex[var_i][var_j][var_k], Ey[var_i][var_j][var_k], Ez[var_i][var_j][var_k]);
					}
					puts("");
				}
				puts("");
			}
			puts("");
			
			// 電界の保存
			saving_electric_field();

			// 磁界の計算
			calc_hfield();

			MPI_Sendrecv( &Hy[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
				&Hy[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);
			MPI_Sendrecv( &Hz[xmax-1][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, right, tag_send,
				&Hz[0][0][0], (ymax+1)*(zmax+1), MPI_DOUBLE, left, tag_recv, MPI_COMM_WORLD, &status);

			// フィールドの出力
			output_field (dir_name[dir_count]);

			// ポインティングパワー計算と出力
#if _EXITATION_FUNC
			calc_poynting_power();
#else
			calc_poynting_powerHz();
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
			//etime == time;//☆計算終了時刻をホールド
			//fprintf(fpparameter, "End Time:\t %s\n", time); 	//計算終了時刻の出力
			//fprintf(fpparameter, "Calculation Time:\t %s\n", time-stime); 	//☆計算時間の出力

			//時刻の出力
		//	e_time = MPI_Wtime();
			//printf ("\ntime = %f\n", e_time - s_time);
			//fprintf(fpparameter, "Calculation Time:\t %s\n", time);
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

	//baba lab
	//_mkdir(strcpy(dir_name, dir_name_def)); 		// ?U?蕪????????e?X?g

	//kuramitsu lab
	mkdir(strcpy(dir_name, dir_name_def), 0755); 		// ?U?蕪????????e?X?g


	if (irank == IRANK_MIN){
		fpparameter = fopen (strcat(strcpy(dir_name, dir_name_def), "/Parameter.txt"), "w");
		allmodel_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_xy.txt"), "w");
		fpallepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/All_Epsilon_xy.txt"), "w");
	}

	if(irank == intObseInPortNum){ // 入射
		allmodel_yz1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_yz1.txt"), "w");
	}
	if(irank == intObseOutPortNum){ // 出射
		allmodel_yz4 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_yz4.txt"), "w");
	}
	if(irank == intObseCenPortNum){			// 中央
		allmodel_yz7 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AllModel_yz7.txt"), "w");
	}

	model_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy.txt"), "w"); 		// 振り分けできるかテスト
	model_xy2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xy2.txt"), "w");
	model_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_yz.txt"), "w");
	model_xz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Model_xz.txt"), "w");

	//fpfez_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Ez_xy.txt"), "w");
	//fpfez2_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Ez2_xy.txt"), "w");
	//fpfez_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Ez_yz.txt"), "w");
	//fpfez_zx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Ez_zx.txt"), "w");
	//fpfhz_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Hz_xy.txt"), "w");
	//fpfhz2_xy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Hz2_xy.txt"), "w");
	//fpfhz_yz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Hz_yz.txt"), "w");
	//fpfhz_zx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Field_Hz_zx.txt"), "w");

	fpepsilonx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_xy.txt"), "w");
	fpepsilony = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz.txt"), "w");
	fpepsilonz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx.txt"), "w");
	fpepsilony2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_yz2.txt"), "w");
	fpepsilonz2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Epsilon_zx2.txt"), "w");

	//fpex = fopen (strcat(strcpy(dir_name, dir_name_def), "/Ex.txt"), "w");
	//fphx = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hx.txt"), "w");
	//fpey = fopen (strcat(strcpy(dir_name, dir_name_def), "/Ey.txt"), "w");
	//fphy = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hy.txt"), "w");
	//fpez = fopen (strcat(strcpy(dir_name, dir_name_def), "/Ez.txt"), "w");
	//fphz = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hz.txt"), "w");

	if(irank == intObseInPortNum){
		fppoynt1 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower1.txt"), "w");
		avpoynt1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePoyntingPower1.txt"), "w");
		fpHz1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hz01.txt"), "w");
	}
	//fppoynt2 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower2.txt"), "w");
	//fppoynt3 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower3.txt"), "w");
	//fppoynt4 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower4.txt"), "w");
	if(irank == intObseOutPortNum){
		fppoynt5 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower5.txt"), "w");
		avpoynt5 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePoyntingPower5.txt"), "w");
		fpHz5 = fopen (strcat(strcpy(dir_name, dir_name_def), "/Hz05.txt"), "w");
	}
	//fppoynt6 = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower6.txt"), "w"); //ポインティン
	//fppoynt1h = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower1h.txt"), "w");
	//fppoynt2h = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower2h.txt"), "w");
	//fppoynt3h = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower3h.txt"), "w");
	//fppoynt4h = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower4h.txt"), "w");
	//fppoynt5h = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower5h.txt"), "w");
	//fppoynt6h = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPower6h.txt"), "w");
	//fppoynt_para = fopen(strcat(strcpy(dir_name, dir_name_def), "/PoyntingPara.txt"), "w"); 		//ポインティングパラメータ保存ファイル
	//fppowerHz1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz1.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz2 = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz2.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz3 = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz3.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz4 = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz4.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz5 = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz5.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz6 = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz6.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz1h = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz1h.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz2h = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz2h.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz3h = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz3h.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz4h = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz4h.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz5h = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz5h.txt"), "w"); 		//Hz保存ファイル
	//fppowerHz6h = fopen (strcat(strcpy(dir_name, dir_name_def), "/PowerHz6h.txt"), "w"); 		//Hz保存ファイル

	//avhz1 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePowerHz1.txt"), "w");
	//avhz5 = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePowerHz5.txt"), "w");

	//avpoynt1h = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePoyntingPower1h.txt"), "w");
	//avpoynt5h = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePoyntingPower5h.txt"), "w");
	//avhz1h = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePowerHz1h.txt"), "w");
	//avhz5h = fopen (strcat(strcpy(dir_name, dir_name_def), "/AveragePowerHz5h.txt"), "w");
}


/*出力用ファイルを閉じる*/
void file_close(){

	if (irank == IRANK_MIN){
		fclose(fpparameter);
		//fclose(allmodel_xy);
		fclose(fpallepsilonx);
	}

	fclose(model_xy);
	fclose(model_xy2);
	fclose(model_yz);
	fclose(model_xz);

	fclose(fpepsilonx);
	fclose(fpepsilony);
	fclose(fpepsilonz);
	fclose(fpepsilony2);
	fclose(fpepsilonz2);

	if(irank == intObseInPortNum){
		fclose(fppoynt1);
		fclose(avpoynt1);
		fclose(fpHz1);
	}
	if(irank == intObseOutPortNum){
		fclose(fppoynt5);
		fclose(avpoynt5);
		fclose(fpHz5);
	}
}


// 計算用パラメータの設定と出力
void parameter(char* dir_name){

	if (irank == IRANK_MIN){

		//printf("Nodes = %d\n", NODE);
		//printf("Cell Size[nm] = %d\n", CELL_SIZE);
		//printf("Time Step[s] = %e\n", dt);
		//printf("Final Time Step: %d\n", Nmax);

		//以下parameter.txtの中身
		fprintf(fpparameter, "XMAX_ALL = %d\n", XMAX_ALL);
		fprintf(fpparameter, "YMAX_ALL = %d\n", YMAX_ALL);
		fprintf(fpparameter, "ZMAX_ALL = %d\n", ZMAX_ALL);
		fprintf(fpparameter, "BOUNDARYLINE = %d\n", BOUNDARYLINE);//15/12/24宮坂
		fprintf(fpparameter, "Nodes = %d\n", NODE);
		fprintf(fpparameter, "Cell Size [nm] = %d\n", CELL_SIZE);
		fprintf(fpparameter, "Time Step [s] = %e\n", dt);
		fprintf(fpparameter, "Final Time Step = %d\n", Nmax);
		fprintf(fpparameter, "Final Time [s] = %e\n", (double) Nmax * dt);
		//☆☆パワーの平均の算出を開始する時間ステップ  (最終計算ステップからの差)
		fprintf(fpparameter, "The starting time of poyinting vector =%e\n",(double) (Nmax- Tcut)*dt);

		//printf("Slab Height [nm] = %d\n", SLAB_HEIGHT);
		//printf("Hole Pitch [nm] = %d\n", PITCH);
		//printf("Hole Diameter [nm] = %d\n", RADIUS*2);
		//printf("Normal PCW Period = %d\n", NORM_PCW_PER);
		//printf("Chirped LSPCW Period = %d\n", CHIRP_3RD_LS_PER);
		//printf("LSPCW Period = %d\n", LSPCW_PER);
		//printf("Hole Column[y] = %d\n", intPcwWid);
		//printf("Hole Row[x] = %d \n", intPcwPer);
		//printf("Hole Start Coordinate = (%d, %d)\n", intPcwStartX, intPcwStartY);
		//dir_name
		fprintf(fpparameter, "Wavelength [nm] = %lf\n", dir_name[1]);
		fprintf(fpparameter, "Upper Clad Index = %lf\n", n_clad);
		fprintf(fpparameter, "Slab Index = %lf\n", n_core);
		fprintf(fpparameter, "Upper Height [nm] = %d\n", CLAD_HEIGHT1);
		fprintf(fpparameter, "Slab Height [nm] = %d\n", SLAB_HEIGHT);
		fprintf(fpparameter, "Hole Pitch [nm] = %d\n", PITCH);
		fprintf(fpparameter, "Pitch Shift Max [nm] = %d\n", PITCH_SHIFT_MAX);
		fprintf(fpparameter, "Hole Diameter [nm] = %d\n", RADIUS*2);
		fprintf(fpparameter, "%d-row Shift [nm] = %d\n", 1, SX1);
		fprintf(fpparameter, "%d-row Shift [nm] = %d\n", 2, SX2);
		fprintf(fpparameter, "%d-row Shift [nm] = %d\n", 3, SX3);
		fprintf(fpparameter, "%Y-direction All Shift [nm] = %d\n", SY);
		//fprintf(fpparameter, "%d-row Shift [nm] = %d\n", LSPCW_ROW, SX3);
		fprintf(fpparameter, "Normal PCW Period = %d\n", NORM_PCW_PER);
		fprintf(fpparameter, "Chirped LSPCW Period = %d\n", CHIRP_3RD_LS_PER);
		fprintf(fpparameter, "Pitch Shift PCW Period = %d\n", PITCH_SHIFT_PER);//☆
		fprintf(fpparameter, "Pitch Shift Chirp PCW Period = %d\n", PITCH_SHIFT_CHIRP_PER);
		fprintf(fpparameter, "LSPCW Period = %d\n", LSPCW_PER);
		fprintf(fpparameter, "Hole Column[y] = %d \n", intPcwWid);
		fprintf(fpparameter, "Hole Row[x] = %d \n", intPcwPer);
		fprintf(fpparameter, "Hole Start Coordinate = (%d, %d)\n", intPcwStartX, intPcwStartY);
		fprintf(fpparameter, "\n");

		fprintf(fpparameter, "Exctation = %d nm (%d cell)\n", EXCT_LEN, EXCT_LEN/CELL_SIZE);
		fprintf(fpparameter, "Exctation to Observation = %d nm (%d cell)\n", EXCT_OBSE_LEN, EXCT_OBSE_LEN/CELL_SIZE);
		fprintf(fpparameter, "Observation to Wire = %d nm (%d cell)\n", OBSE_WIRE_LEN, OBSE_WIRE_LEN/CELL_SIZE);
		fprintf(fpparameter, "Output Wire to Termination = %d nm (%d cell)\n", WIRE_OUTPUT_LEN, WIRE_OUTPUT_LEN/CELL_SIZE);
		fprintf(fpparameter, "Termination Length = nm %d (%d cell)\n", WIRE_OUTPUT_OFFSET, WIRE_OUTPUT_OFFSET/CELL_SIZE);
		if(intWireWid < 0){
			fprintf(fpparameter, "Wire Width = %d nm (%d cell)\n", (intWireWid + 1) * CELL_SIZE, ((intWireWid + 1) * CELL_SIZE)/CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "Wire Width = %d nm (%d cell)\n", intWireWid * CELL_SIZE, ((intWireWid + 1) * CELL_SIZE) / CELL_SIZE);
		}
		fprintf(fpparameter, "PCW Slab Termination Length = %d nm (%d cell)\n", PCW_SiSLAB_TERMINATION_LEN, PCW_SiSLAB_TERMINATION_LEN/CELL_SIZE);
		if(PCW_SiSLAB_OFFSET < 0){
			fprintf(fpparameter, "PCW Slab Offset = %d nm (%d cell)\n", PCW_SiSLAB_OFFSET + CELL_SIZE, (PCW_SiSLAB_OFFSET + CELL_SIZE)/CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "PCW Slab Offset = %d nm (%d cell)\n", PCW_SiSLAB_OFFSET, PCW_SiSLAB_OFFSET/CELL_SIZE);
		}
		if(PCW_WIDTH_CHIRP < 0){//入出射 //負になることは基本ないはずであるが？
			fprintf(fpparameter, "PCW Width Chirp_IN, OUT = %d nm (%d cell), %d nm (%d cell)\n", PCW_WIDTH_CHIRP + CELL_SIZE, (PCW_WIDTH_CHIRP + CELL_SIZE)/CELL_SIZE, PCW_WIDTH_CHIRP_OUT+CELL_SIZE, (PCW_WIDTH_CHIRP_OUT + CELL_SIZE)/CELL_SIZE);
		}
		else{
			fprintf(fpparameter, "PCW Width Chirp_IN, OUT = %d nm (%d cell), %d nm (%d cell)\n", PCW_WIDTH_CHIRP, PCW_WIDTH_CHIRP/CELL_SIZE, PCW_WIDTH_CHIRP_OUT, PCW_WIDTH_CHIRP_OUT/CELL_SIZE);
		}

		/*if (PCW_WIDTH_CHIRP_OUT < 0) {//出射//入出射まとめたので要らない
			fprintf(fpparameter, "PCW Width Chirp_OUT = %d nm (%d cell),", PCW_WIDTH_CHIRP_OUT + CELL_SIZE, (PCW_WIDTH_CHIRP_OUT + CELL_SIZE)/CELL_SIZE);
		}
		else {
			fprintf(fpparameter, "PCW Width Chirp_OUT = %d nm (%d cell)", PCW_WIDTH_CHIRP_OUT, PCW_WIDTH_CHIRP_OUT/CELL_SIZE);
		}*/

		/*if (PCW_WIDTH_CHIRP == PCW_WIDTH_CHIRP_OUT){
			fprintf(fpparameter, " symmetric Wavelength Width Chirp\n");
		}
		else{
			fprintf(fpparameter, " symmetric Wavelength Width Chirp\n");
		}*/


		if(LSPCW_SHIFT_DESCRETE == FALSE){
			fprintf(fpparameter, "LSPCW_SHIFT_DESCRETE = FALSE\n");
		}
		else{
			fprintf(fpparameter, "LSPCW_SHIFT_DESCRETE = TRUE\n");
		}

		fprintf(fpparameter, "y = %d (light mode profile(zx dimenion)) \n", y_zx);

		fprintf(fpparameter, "\n");
	}

	// 励振関数定数の設定
	lambda = atof(dir_name) * 1e-9;
	omega0 = 2.0*PI*C0/lambda;
	sigma = omega0 * delta_omega;
} //420

/*配列の初期化*/
void initialize_matrix(){

	//各ノードの座標
	if(irank != IRANK_MAX){
		xmax = XMAX;
		ymax = YMAX;
		zmax = ZMAX;
	}

	//最後のノードだけのりしろ不要なのでx方向に1セル小さい
	if(irank == IRANK_MAX){
		xmax = XMAX - 1;
		ymax = YMAX;
		zmax = ZMAX;
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
			for(z = 0; z < zmax+1; z++){
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
			for(z = 0; z < zmax+1; z++){
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
			for(z = 0; z < zmax; z++){
				// セルの目印
				cell[x][y][z] = 0;
			}
		}
	}


	/****************************** Murの吸収境界条件 ******************************/

	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax+1; z++){
			Exn2y00[x][z] = 0.0;
			Exn1y00[x][z] = 0.0;
			Exn2y01[x][z] = 0.0;
			Exn1y01[x][z] = 0.0;
			//exn2ym1[x][z] = 0.0;
			//exn1ym1[x][z] = 0.0;
			//exn2ym0[x][z] = 0.0;
			//exn1ym0[x][z] = 0.0;
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax+1; y++){
			Exn2z00[x][y] = 0.0;
			Exn1z00[x][y] = 0.0;
			Exn2z01[x][y] = 0.0;
			Exn1z01[x][y] = 0.0;
			//exn2zm1[x][y] = 0.0;
			//exn1zm1[x][y] = 0.0;
			//exn2zm0[x][y] = 0.0;
			//exn1zm0[x][y] = 0.0;
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Eyn2z00[x][y] = 0.0;
			Eyn1z00[x][y] = 0.0;
			Eyn2z01[x][y] = 0.0;
			Eyn1z01[x][y] = 0.0;
			//eyn2zm1[x][y] = 0.0;
			//eyn1zm1[x][y] = 0.0;
			//eyn2zm0[x][y] = 0.0;
			//eyn1zm0[x][y] = 0.0;
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax+1; z++){
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
		for(z = 0; z < zmax; z++){
			Ezn2y00[x][z] = 0.0;
			Ezn1y00[x][z] = 0.0;
			Ezn2y01[x][z] = 0.0;
			Ezn1y01[x][z] = 0.0;
			//ezn2ym1[x][z] = 0.0;
			//ezn1ym1[x][z] = 0.0;
			//ezn2ym0[x][z] = 0.0;
			//ezn1ym0[x][z] = 0.0;
		}
	}
	for(y = 0; y <= ymax; y++){
		for(z = 0; z <= zmax-1; z++){
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

				//n_temp = CLAD;
				//epsilon_temp = EPSILON0;
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

		//struct PNUM Pnum[intPcwWid*10][intPcwPer*10]; 	// 円柱の中心座標
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

				//if(y % 2 == 1){		// if y:even
				//	Pnum[0][y].X = intWireLen1 + intPcwStartX - 1;	// 配列の引数に使用するので-1
				//}
				//else{				// if y:odd 0.5Aずらす
				//	Pnum[0][y].X = intWireLen1 + intPcwStartX + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
				//}

				//if(y != 0){
				//	Pnum[0][y].Y = Pnum[0][y-1].Y + intPitchY; 		//ｙ座標を(root3)/2*intPitchXだけずらす
				//}
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
				input_PitchShiftPcw_Xend = input_NormPcw_Xend + intPitchShiftPcwPer + intPitchShiftChirpPcwPer; //Nに相当
				input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend; //Nに相当
				//input_PitchShiftPcw_Xend = input_NormPcw_Xend + intPitchShiftPcwPer;
				//if (intPitchShiftChirpPcwPer > 0){
				//	input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend + (intPitchShiftChirpPcwPer - 1);
				//}
				//else{
				//	input_PitchShiftChirpPcw_Xend = input_PitchShiftPcw_Xend;
				//}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/

				if (intChirp3rdLsPer > 0){
					input_Chirp_Ls_Xend = input_PitchShiftChirpPcw_Xend + (intChirp3rdLsPer);
				}
				else{
					input_Chirp_Ls_Xend = input_PitchShiftChirpPcw_Xend;
				}
				Lspcw_Xend = input_Chirp_Ls_Xend + intLspcwPer + 1; //テーパ周期+LSPCW周期+1
				if (intChirp3rdLsPer > 0){
					output_Chirp_Ls_Xend = Lspcw_Xend + (intChirp3rdLsPer - 1);
				}
				else{
					output_Chirp_Ls_Xend = Lspcw_Xend;
				}
				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend; //テーパ周期+LSPCW周期
				output_PitchShiftPcw_Xend = output_PitchShiftChirpPcw_Xend + intPitchShiftPcwPerOut + intPitchShiftChirpPcwPer;
				//if (intPitchShiftChirpPcwPer > 0){
				//	output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend + (intPitchShiftChirpPcwPer - 1);
				//}
				//else{
				//	output_PitchShiftChirpPcw_Xend = output_Chirp_Ls_Xend;
				//}
				//output_PitchShiftPcw_Xend = output_PitchShiftChirpPcw_Xend + intPitchShiftPcwPer;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/

				output_PCW_Xend = output_Chirp_Ls_Xend + intNormPcwPer;

				int intPitchShiftX = INT_DIV (PITCH_SHIFT_MAX, CELL_SIZE);
				int intPitchShiftY = (INT_DIV((PITCH_SHIFT_MAX * sqrt(3.0)/2 + 0.5), CELL_SIZE));
				int intPCWwidthChirp = INT_DIV (PCW_WIDTH_CHIRP, CELL_SIZE);
				int intPCWwidthChirpOut = INT_DIV (PCW_WIDTH_CHIRP_OUT, CELL_SIZE);
				int intPreviousPCWwidthOffset;
				int intNowPCWwidthOffset;
				int intNowPCWwidthOffset2;

				/******************** 格子シフト量チャープ(2013/7/19) ********************/
				double dblPitchShiftChirpX, dblPitchShiftChirpY, dblPitchShiftChirpY2;
				int intPitchShiftChirpX, intPitchShiftChirpX2, intPitchShiftChirpY;
				/******************** 格子シフト量チャープ(2013/7/19) ********************/




				/****************************** LSPCW ******************************/
						for (x = 0; x < intPcwPer; x++) {//{}の終わりは1740付近 円孔行1, 2, 3…はx = 0, 1, 2…に対応


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
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);//☆ここも変わらない
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
								for (y2 = intPcwWid - 1; y2 >= 0; y2--) {//☆ここも変わらない
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
									//Pnum[x][y2].X -= INT_DIV( (SX1 + (SX1 * (double)(intPitchShiftX / intPitchX) * (x - intPitchShiftChirpPcwPer) / (double) intPitchShiftChirpPcwPer)), CELL_SIZE)-3;
									Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								/******************** 導波路1列目シフト構造(2013/7/12) ********************/

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE);
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;

								///これをいれるとおかしくなる
								/*
							if (y2 == intPcwWid - 2) {
									Pnum[x][y2].Y = Pnum[x - 1][intPcwWid - 1].Y - (int)(((double)intPitchShiftY - dblPitchShiftChirpY * (x - intPitchShiftPcwPer)) * (intPcwWid - 1 - y2)) + 1;//	☆ここ動かすと変わる　y方向に動かしたいセル数
								}
								*/
							}
						}
						else{
							/******************** 格子シフト量チャープ(2013/7/19) ********************/
							intPitchShiftChirpX = intPitchShiftX; //19
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
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchShiftY - dblPitchShiftChirpY * (x - intPitchShiftPcwPer)) * (intPcwWid-1 - y2));//


										if (y2 == intPcwWid - 2) {//☆ -のあとに動かしたいPCWの列(中心が1) 1だと動作しない

											if(x == 3){//☆ ==のあとに動かしたいPCWの行(手前が0)．ただ0の場合はここではできない
												Pnum[x][y2].Y = Pnum[x - 1][intPcwWid - 1].Y - (int)(((double)intPitchShiftY - dblPitchShiftChirpY * (x - intPitchShiftPcwPer)) * (intPcwWid - 1 - y2))  ;//	☆ここ動かすと変わる　y方向に動かしたいセル数
											}
										}



									}
									else{
										Pnum[x][y2].Y = Pnum[x-1][y2].Y;//☆ここ動かすと変わる
									}
								}
								else{
									if (y2 < intPcwWid-1){
										Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);//☆
									}
									else{
										Pnum[x][y2].Y = Pnum[x - 1][y2].Y; //☆
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
								//s_x3 = INT_DIV( (SX3 + (SX3 * (double)(intPitchShiftX / intPitchX) * (x - intPitchShiftChirpPcwPer) / (double) intPitchShiftChirpPcwPer)), CELL_SIZE);
								//s_x3 = INT_DIV( (SX3 * (double)(intPitchShiftX / intPitchX)), CELL_SIZE);
								/******************** 格子シフト量チャープ(2013/7/19) ********************/

							}
							if (y == intPcwWid - 4){
								s_x4 = INT_DIV(SX4, CELL_SIZE);
							}
						}
					}

					// 入射 チャープLSPCW ☆☆1/27 入出射の円孔を少し動かしてみる
					else if (x >= input_NormPcw_Xend && x < input_Chirp_Ls_Xend){
						if (x == 0){
							y_poo = 0; y_poo2 = 0;
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchY * y_poo; //☆		//ｙ座標を(root3)/2*intPitchXだけずらす



								if(y2 % 2 == 1){		// if y:even
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X + intPitchX * x - 1;	// 配列の引数に使用するので-1
								}
								else{				// if y:odd 0.5Aずらす
									Pnum[x][y2].X = Pnum_Init[0][intPcwWid-1].X+ intPitchX * x + INT_DIV (intPitchX, 2.0) - 1;	// 配列の引数に使用するので-1
								}

								/******************** 導波路列全体シフト(2013/7/12) ********************/
								Pnum[x][y2].Y -= INT_DIV(SY, CELL_SIZE)+100; //☆ここも変わらない
								/******************** 導波路列全体シフト(2013/7/12) ********************/

								y_poo++;
							}
						}
						else{
							for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y = Pnum[x-1][y2].Y;//ここをいじっても変わらない．
								Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchX;//

								/*
								if (y2 == 1) {//☆
									Pnum[x][y2].Y = Pnum[x - 1][y2].Y +10;
								}*/
							}

							if (y == intPcwWid - 3){
								if (intChirp3rdLsPer == 0){
									s_x3 = 0;//☆これをいじっても×
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
								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								//Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchShiftY * y_poo; 		//ｙ座標を(root3)/2*intPitchXだけずらす

								//if(y2 % 2 == 1){		// if y:even
								//	Pnum[x][y2].X = Pnum[x-1][intPcwWid-1].X + intPitchShiftX;	// 配列の引数に使用するので-1
								//}
								//else{				// if y:odd 0.5Aずらす
								//	Pnum[x][y2].X = Pnum[x-1][intPcwWid-1].X + intPitchShiftX + INT_DIV (intPitchX, 2.0);	// 配列の引数に使用するので-1
								//}
								/******************** 格子シフト量チャープ(2013/7/19) ********************/


								/******************** 導波路1列目シフト構造(2013/7/12) ********************/
								if (y2 != intPcwWid - 1){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
									//Pnum[x][y2].X -= INT_DIV( (SX1 + (SX1 * (double)(intPitchShiftX / intPitchX) * (x - output_PitchShiftChirpPcw_Xend) / (double) intPitchShiftChirpPcwPer)), CELL_SIZE);
									//Pnum[x][y2].X -= INT_DIV(SX1, CELL_SIZE);
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

							//if (x > output_PitchShiftChirpPcw_Xend){
							//	dblPitchShiftChirpY = (intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
							//}
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

							//if(x > output_PitchShiftChirpPcw_Xend){
							//if (x >= output_PitchShiftChirpPcw_Xend){
							//	if (y2 < intPcwWid-1){
							//		Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - (int) ( ((double)intPitchY + dblPitchShiftChirpY * (x - output_PitchShiftChirpPcw_Xend + 1)) * (intPcwWid-1 - y2));
							//	}
							//	else{
							//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;
							//	}
							//}
							//else{
							//	if (y2 < intPcwWid-1){
							//		Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);
							//	}
							//	else{
							//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;
							//	}
							//}
							//}
							//else{
							//	if (y2 < intPcwWid-1){
							//		Pnum[x][y2].Y = Pnum[x-1][intPcwWid-1].Y - intPitchY * (intPcwWid-1 - y2);
							//	}
							//	else{
							//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;
							//	}
							//}
							Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						}

						////else{
						//intPitchShiftChirpX = intPitchX;
						//intPitchShiftChirpX2 = intPitchX;
						//intPitchShiftChirpY = 0;
						//if (x >= output_PitchShiftChirpPcw_Xend){
						//	intPitchShiftChirpX += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftChirpPcwPer)));
						//	intPitchShiftChirpX2 += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftChirpPcwPer - 1)));
						//	//if (x > output_PitchShiftChirpPcw_Xend){
						//	dblPitchShiftChirpY = -(intPitchShiftY - intPitchY) / (double) intPitchShiftChirpPcwPer;
						//	dblPitchShiftChirpY2 = -(intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 1);
						//	intPitchShiftChirpY = (int) ((intPitchShiftY - intPitchY) * (output_PitchShiftChirpPcw_Xend-x-1) / (double) intPitchShiftChirpPcwPer) - (int) ((intPitchShiftY - intPitchY) * (output_PitchShiftChirpPcw_Xend-x) / (double) intPitchShiftChirpPcwPer);
						//	//}
						//}

						//for (y2 = intPcwWid-1; y2 >= 0; y2--){
						//	if (y2 < intPcwWid-1){
						//		if (y2 % 2 == 1){
						//			Pnum[x][y2].Y = Pnum[x-1][y2].Y + (int) (dblPitchShiftChirpY * (double) (intPcwWid-1 - y2));
						//		}
						//		else{
						//			Pnum[x][y2].Y = Pnum[x-1][y2].Y + intPitchShiftChirpY + (int) (dblPitchShiftChirpY2 * (double) (intPcwWid-1 - y2 - 1));
						//		}
						//	}
						//	else{
						//		Pnum[x][y2].Y = Pnum[x-1][y2].Y;
						//	}
						//	if (y2 % 2 == 1){
						//		Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						//	}
						//	else{
						//		Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX2;
						//	}
						//}

						////for (y2 = intPcwWid-1; y2 >= 0; y2--){
						////	Pnum[x][y2].Y = Pnum[x-1][y2].Y;
						////	Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftX;
						////}
						////}
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
								//s_x3 = INT_DIV( (SX3 + (SX3 * (double)(intPitchShiftX / intPitchX) * (x - output_PitchShiftChirpPcw_Xend) / (double) intPitchShiftChirpPcwPer)) , CELL_SIZE);
								//s_x3 = INT_DIV( (SX3 * (double)(intPitchShiftX / intPitchX)), CELL_SIZE);
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
						///******************** 格子シフト量チャープ(2013/7/19) ********************/
						////// 格子定数変化PCWとの出射接続部
						//if (x == output_Chirp_Ls_Xend - 1){
						//	intPitchShiftChirpX = intPitchX;
						//	intPitchShiftChirpY = 0;

						//	intPitchShiftChirpX += (int) ((intPitchShiftX - intPitchX) * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftChirpPcwPer - 1)));
						//	//if (x > output_PitchShiftChirpPcw_Xend){
						//	dblPitchShiftChirpY = -(intPitchShiftY - intPitchY) / (double) intPitchShiftChirpPcwPer;
						//	dblPitchShiftChirpY2 = -(intPitchShiftY - intPitchY) / (double) (intPitchShiftChirpPcwPer - 2);
						//	intPitchShiftChirpY = (int) ((intPitchShiftY - intPitchY) * (intPitchShiftChirpPcwPer-x-1) / (double) intPitchShiftChirpPcwPer) - (int) ((intPitchShiftY - intPitchY) * (intPitchShiftChirpPcwPer-x) / (double) intPitchShiftChirpPcwPer);
						//	//}

						//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
						//		if (y2 < intPcwWid-1){
						//			if (y2 % 2 == 1){
						//				//Pnum[x][y2].Y = Pnum[x-1][y2].Y + (int) (dblPitchShiftChirpY * (double) (intPcwWid-1 - y2));
						//			}
						//			else{
						//				Pnum[x][y2].Y = Pnum[x-1][y2].Y + intPitchShiftChirpY + (int) (dblPitchShiftChirpY2 * (double) (intPcwWid-1 - y2 - 1.5));
						//			}
						//		}
						//		else{
						//			Pnum[x][y2].Y = Pnum[x-1][y2].Y;
						//		}
						//		Pnum[x][y2].X = Pnum[x-1][y2].X + intPitchShiftChirpX;
						//	}
						//	//	y_poo = 0;
						//	//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
						//	//		if(y2 % 2 == 0){
						//	//			Pnum[x][y2].X += (intPitchShiftX - intPitchX);
						//	//			Pnum[x][y2].Y -= (intPitchShiftY - intPitchY) * y_poo;
						//	//		}
						//	//		y_poo++;
						//	//	}
						//}
						///******************** 格子シフト量チャープ(2013/7/19) ********************/


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
					Pnum[x][y2].X += 1;
					if (PCW_WIDTH_CHIRP != 0){
						if (x < input_PitchShiftPcw_Xend){ //テーパ周期数はNより1小さい

							/////格子定数チャープ入射ここから(田村)
							if(PITCH_SHIFT_CHIRP_PER2 > 0 && x >= 0){
								if(PITCH_SHIFT_CHIRP_PER2 >= x){
									y_poo = 0;
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										if(y2 % 2 == 1 && x != 0){
											Pnum[x][y2].X += INT_DIV((PITCH_SHIFT_MAX2 - PITCH)/PITCH_SHIFT_CHIRP_PER2*(PITCH_SHIFT_CHIRP_PER2 - x + 2),CELL_SIZE );
										}
										if(y2 % 2 == 0 && x != 0){
											Pnum[x][y2].X += INT_DIV((PITCH_SHIFT_MAX2 - PITCH)/PITCH_SHIFT_CHIRP_PER2*(PITCH_SHIFT_CHIRP_PER2 - x + 2),CELL_SIZE );
										}
										if(y2 != intPcwWid-1){
											Pnum[x][y2].Y -= INT_DIV(((PITCH_SHIFT_MAX2 * sqrt(3.0) / 2) - (PITCH * sqrt(3.0) / 2))/PITCH_SHIFT_CHIRP_PER2*(PITCH_SHIFT_CHIRP_PER2 - x )* y_poo,CELL_SIZE) ;
										}
										y_poo++;
									}
								}
							}
							/////格子定数チャープ入射ここまで

							////////入射シフト量チャープここから(田村)
							if(intChirp2ndLsPer != 0){
							//3列目
							if (y == intPcwWid - 3){
								if (intChirp2ndLsPer == 0){
									s_x3 = 0;
								}
								else if (INT_DIV ((SX3 / intChirp2ndLsPer) * (x +5 - input_PitchShiftPcw_Xend),CELL_SIZE)<INT_DIV(SX3,CELL_SIZE)){
									s_x3 = INT_DIV ((SX3 / intChirp2ndLsPer) * (x +5 - input_PitchShiftPcw_Xend),CELL_SIZE);
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
								else if (INT_DIV ((SX2 / intChirp2ndLsPer) * (x +5 - input_PitchShiftPcw_Xend),CELL_SIZE)<INT_DIV(SX2,CELL_SIZE)){
									s_x2 = INT_DIV ((SX2 / intChirp2ndLsPer) * (x +5 - input_PitchShiftPcw_Xend),CELL_SIZE);
								}
								else{
									s_x2 = INT_DIV(SX2,CELL_SIZE);
								}
							}
							}
							////////入射シフト量チャープここまで
							if (x == 0){ //円孔1行目
								intNowPCWwidthOffset = intPCWwidthChirp; //wtaperのセル数
								for (y2 = intPcwWid-1; y2 >= 0; y2--){
									Pnum[x][y2].Y -= intNowPCWwidthOffset;
								}
								intPreviousPCWwidthOffset = intNowPCWwidthOffset; //wtaperのセル数
							}
							else{ //円孔2行目以降
								// PCW_WIDTH_CHIRP用微調整☆☆ 11/9 ここがおそらく幅チャ
								intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * (x / (double) (input_PitchShiftPcw_Xend-1)));// 右辺 = wtaperのセル数(1 - x/N-1)
								intNowPCWwidthOffset2 = intPCWwidthChirp - (int)(intPCWwidthChirp * (x / (double)(input_PitchShiftPcw_Xend -2.5)));//下式によりテーパの分割が4+1になる //N =9で2.5

								//intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * ((x + 1) / (double) (input_PitchShiftPcw_Xend-1)));
								if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){ //intPreviousPCWwidthOffsetは多分上のif文のintNowPCWwidthOffset(wtaperのセル数)に相当
									//→x = 0以外ではif文を通過できるはず
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
								Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset); // (1/6)wtaperずつ //オリジナル

									//// 2016/11/9作成
									//	if (x >= input_PitchShiftPcw_Xend -4) { //x＞N-4で幅チャープを適用
									//		Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset2);
									//	}else{
									//		Pnum[x][y2].Y = Pnum[x][y2].Y;
									//	}//

									}
									intPreviousPCWwidthOffset = intNowPCWwidthOffset;
								}
							}
						}
						else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){

							/////格子定数チャープ出射ここから(田村氏)
							if(PITCH_SHIFT_CHIRP_PER2_OUT > 0 && x >= 0){
								if(PITCH_SHIFT_CHIRP_PER2_OUT + output_PitchShiftChirpPcw_Xend+1 >= x){
									y_poo = 0;
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										if(x != 0){
											Pnum[x][y2].X += INT_DIV((PITCH_SHIFT_MAX2 - PITCH)/PITCH_SHIFT_CHIRP_PER2_OUT*(x - output_PitchShiftChirpPcw_Xend + 1),CELL_SIZE );
										}
										if(y2 != intPcwWid-1 ){
											//if(y2%2 == 0)
											//Pnum[x][y2].Y -= INT_DIV(((PITCH_SHIFT_MAX2 * sqrt(3.0)/2) - (PITCH * sqrt(3.0)/2))/PITCH_SHIFT_CHIRP_PER2*(x - output_PitchShiftChirpPcw_Xend +1)* y_poo,CELL_SIZE) ;
											//if(y2%2 == 1)
											Pnum[x][y2].Y -= INT_DIV(((PITCH_SHIFT_MAX2 * sqrt(3.0)/2) - (PITCH * sqrt(3.0)/2))/PITCH_SHIFT_CHIRP_PER2_OUT*(x - output_PitchShiftChirpPcw_Xend +1)* y_poo,CELL_SIZE);
											y_poo++;
										}
										//y_poo++;
									}
								}
							}


							/////格子定数チャープ出る射ここまで

							////////出射シフト量チャープここから(田村)
							if(intChirp2ndLsPer != 0){
							//3列目
							if (y == intPcwWid - 3){
								if (intChirp2ndLsPer == 0){
									s_x3 = 0;
								}
								else if(INT_DIV ((SX3 / intChirp2ndLsPerOut) * (output_PitchShiftPcw_Xend - x-2),CELL_SIZE)<INT_DIV(SX3,CELL_SIZE)){
									s_x3 = INT_DIV ((SX3 / intChirp2ndLsPerOut) * (output_PitchShiftPcw_Xend - x-2),CELL_SIZE);
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
								else if( INT_DIV ((SX2 / intChirp2ndLsPerOut) * (output_PitchShiftPcw_Xend - x-2),CELL_SIZE)<INT_DIV(SX2,CELL_SIZE)){
									s_x2 = INT_DIV ((SX2 / intChirp2ndLsPerOut) * (output_PitchShiftPcw_Xend - x-2),CELL_SIZE);
								}
								else{
									s_x2 = INT_DIV(SX2,CELL_SIZE);
								}
							}
							}
						//////////出射シフト量チャープここまで

							if (x == output_PitchShiftChirpPcw_Xend){ //全LSPCW行数 30～
								//for (y2 = intPcwWid-1; y2 >= 0; y2--){
								//	Pnum[x][y2].Y -= intNowPCWwidthOffset;
								//}
								intPreviousPCWwidthOffset = 0;

								intNowPCWwidthOffset = (int) (intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double) (intPitchShiftPcwPerOut-1)) + 0.9);
								if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
									for (y2 = intPcwWid-1; y2 >= 0; y2--){
										if (y2 % 2 == 0){//これを付けないと出射端の位置がずれる
											Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
										}
									}
									intPreviousPCWwidthOffset = intNowPCWwidthOffset;
								}
								intPreviousPCWwidthOffset = 0;
							}
							else{
								// PCW_WIDTH_CHIRP用微調整
								//intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x - output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)) + 0.8);


								for (y2 = intPcwWid-1; y2 >= 0; y2--){
										if (y2 % 2 == 1) {
											intNowPCWwidthOffset = (int)(intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend) / (double)(intPitchShiftPcwPerOut - 1)) + 0.9);
											intNowPCWwidthOffset2 = (double)(intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend) / (double)(intPitchShiftPcwPerOut - 4.8))+0.1);
										}
										else {
											intNowPCWwidthOffset = (int)(intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double)(intPitchShiftPcwPerOut - 1)) + 0.9);
											intNowPCWwidthOffset2 = (double)(intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend + 1) / (double)(intPitchShiftPcwPerOut - 4.8))+0.1);
										}



									if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
										Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset); //オリジナル

										//// 16/11/10作成
										//if (x<= output_PitchShiftChirpPcw_Xend+4) {
										//	Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset2);
										//	if (y2 % 2 == 0) {
										//		if (x == output_PitchShiftChirpPcw_Xend + 3) {
										//			Pnum[x][y2].Y = Pnum[x][y2].Y+2;//
										//		}
										//	}

										//}
										//else {
										//	if ( y2 % 2 == 0) {
										//		if (x == output_PitchShiftChirpPcw_Xend + 4){
										//			Pnum[x][y2].Y = Pnum[x][y2].Y;//
										//	}
										//	}else{
										//		Pnum[x][y2].Y = Pnum[x][y2].Y;
										//	}
										//}//

										/*if (PCW_WIDTH_Para == 1 && x==27){
												Pnum[x][y2].Y = Pnum[x][y2].Y-1;
										}
										if (PCW_WIDTH_Para == 1 && x==29){
												Pnum[x][y2].Y = Pnum[x][y2].Y+1;
										}
										if (PCW_WIDTH_Para == 1 && x==28){
												Pnum[x][y2].Y = Pnum[x][y2].Y+1;
										}*/

									}
								}
								//intPreviousPCWwidthOffset = (int) (intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend) / (double) (intPitchShiftPcwPerOut-1)) + 0.9);
								intPreviousPCWwidthOffset = (int)(intPCWwidthChirpOut * ((x - output_PitchShiftChirpPcw_Xend) / (double)(intPitchShiftPcwPerOut - 1)) + 0.9);//これは変えない
							}
						}
					}//幅チャ+格チャ終




					//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//if (x < input_PitchShiftPcw_Xend){
					//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
					//		if (x == 0){
					//				intNowPCWwidthOffset = intPCWwidthChirp;
					//				Pnum[x][y2].Y -= intNowPCWwidthOffset;
					//				intPreviousPCWwidthOffset = intPCWwidthChirp;
					//			}
					//		}
					//		else{
					//			if (PCW_WIDTH_CHIRP != 0){
					//				intNowPCWwidthOffset = intPCWwidthChirp - (int) (intPCWwidthChirp * (x / (double) (input_PitchShiftPcw_Xend-1)));
					//				if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
					//					Pnum[x][y2].Y -= (-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
					//				}
					//			}
					//		}
					//	}
					//}
					//else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
					//	for (y2 = intPcwWid-1; y2 >= 0; y2--){
					//		if (x == output_PitchShiftChirpPcw_Xend){
					//			if (PCW_WIDTH_CHIRP != 0){
					//				intNowPCWwidthOffset = intPCWwidthChirp;
					//				Pnum[x][y2].Y -= intNowPCWwidthOffset;
					//				intPreviousPCWwidthOffset = intPCWwidthChirp;
					//			}
					//		}
					//		else{
					//			if (PCW_WIDTH_CHIRP != 0){
					//				intNowPCWwidthOffset = (int) (intPCWwidthChirp * ((x-output_PitchShiftChirpPcw_Xend) / (double) (input_PitchShiftPcw_Xend-1)));
					//				if ( abs(intPreviousPCWwidthOffset) != abs(intNowPCWwidthOffset) ){
					//					Pnum[x][y2].Y -= -(-intPreviousPCWwidthOffset + intNowPCWwidthOffset);
					//				}
					//			}

					//		}
					//	}
					//}



					//}

					// 通常PCWあり
					//else{
					//	if (x != 0){
					//		Pnum[x][y].X = Pnum[x-1][y].X + intPitchX; 		// intPitchXだけ+X座標に置く
					//		Pnum[x][y].Y = Pnum[x-1][y].Y; 					// 同じY座標に置く
					//	}

					//	s_x3 = 0;

					//	if (y == intPcwWid - LSPCW_ROW){

					//		// 入射 通常PCW
					//		if (intNormPcwPer != 0 && x < input_NormPcw_Xend + 1){
					//		}
					//		// 入射 チャープLSPCW
					//		if (x >= input_NormPcw_Xend && x < input_Chirp_Ls_Xend + 1){
					//			if (intChirp3rdLsPer == 0){
					//				s_x3 = 0;
					//			}
					//			else{
					//				s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (x + 1 - input_NormPcw_Xend);
					//			}
					//		}
					//		// LSPCW
					//		if (x >= input_Chirp_Ls_Xend && x < Lspcw_Xend){
					//			s_x3 = intSx3Per;
					//		}
					//		// 出射 チャープLSPCW
					//		if (x >= Lspcw_Xend && x < output_Chirp_Ls_Xend){
					//			if (intChirp3rdLsPer == 0){
					//				s_x3 = 0;
					//			}
					//			else{
					//				s_x3 = INT_DIV (intSx3Per, intChirp3rdLsPer) * (output_Chirp_Ls_Xend - x);
					//			}
					//		}
					//		// 入射 通常PCW
					//		if (x >= output_Chirp_Ls_Xend && x < output_PCW_Xend){
					//		}
					//	}
					//}

					/******************** 格子シフト量チャープ(2013/7/19) ********************/
					if (PITCH_SHIFT_MAX > PITCH && intNormPcwPer == 0 && intPitchShiftPcwPer + intPitchShiftChirpPcwPer != 0){
						//// 入射 格子定数変化PCW
						//if (x < input_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						//// 出射 格子定数変化PCW
						//else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						// 格子定数変化PCWとの出射接続部
						//if (x <= 0 + 1 || x >= intPcwPer - (1+1) ){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}

						//// 入射 格子定数変化PCW
						//if (x < input_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						//// 出射 格子定数変化PCW
						//else if (x >= output_PitchShiftChirpPcw_Xend && x < output_PitchShiftPcw_Xend){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}
						//// 格子定数変化PCWとの出射接続部
						//if (x == output_Chirp_Ls_Xend - 1){
						//	if (y == 0){
						//		continue;	// 解析領域の都合上，格子定数変化ではPCW列数をintPcwWid-1に
						//	}
						//}

						//// 入出射対象構造
						//if (x == output_PitchShiftChirpPcw_Xend - 1) {
						//	if (y % 2 == 1){
						//		continue;
						//	}
						//}
					}
					/******************** 格子シフト量チャープ(2013/7/19) ********************/

					if (SX2 == 0 && SX4 == 0) {
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);

						/*
						{
							if (x <= 6 || x >= 22) {
								mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 12);

								if (x == 22 && y % 2 == 1)
									mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
							}
							else {
								mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
							}
						}*/


						if (x <= 6 || x >= 22) {
						if (x == 0 || x == 29 && y%2 == 1 || x == 28 && y%2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 9);
						if (x == 1 || x == 28 && y % 2 == 1 || x == 27 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 10);
						if (x == 2 || x == 27 && y % 2 == 1 || x == 26 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 11);
						if (x == 3 || x == 26 && y % 2 == 1 || x == 25 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 12);
						if (x == 4 || x == 25 && y % 2 == 1 || x == 24 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 13);
						if (x == 5 || x == 24 && y % 2 == 1 || x == 23 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 14);
						if (x == 6 || x == 23 && y % 2 == 1 || x == 22 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 15);

						if ( x == 22 && y % 2 == 1)
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
						}
						else{
						mcircle(Pnum[x][y].X + s_x3, Pnum[x][y].Y, z, 1);
						}
					}

					if (SX3 == 0 && SX4 == 0) {
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);

						/*
						if (x <= 6 || x >= 22) {
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 12);

						if (x == 22 && y % 2 == 1)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
						}
						else {
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
						}

						/*
						if (x <= 6 || x >= 22) {
						if (x == 0 || x == 29 && y % 2 == 1 || x == 28 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 9);
						if (x == 1 || x == 28 && y % 2 == 1 || x == 27 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 10);
						if (x == 2 || x == 27 && y % 2 == 1 || x == 26 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 11);
						if (x == 3 || x == 26 && y % 2 == 1 || x == 25 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 12);
						if (x == 4 || x == 25 && y % 2 == 1 || x == 24 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 13);
						if (x == 5 || x == 24 && y % 2 == 1 || x == 23 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 14);
						if (x == 6 || x == 23 && y % 2 == 1 || x == 22 && y % 2 == 0)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 15);

						if (x == 22 && y % 2 == 1)
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
						}
						else {
						mcircle(Pnum[x][y].X + s_x2, Pnum[x][y].Y, z, 1);
						}
						*/
					}

					if (SX2 == 0 && SX3 == 0) {
						mcircle(Pnum[x][y].X + s_x4, Pnum[x][y].Y, z, 1);
					}

					// 出射CORE細線導波路との接続部分の調整
					if(x == intPcwPer - 1){
						// 格子定数変化PCWあり
						/******************** 格子シフト量チャープ(2013/7/19) ********************/
						if (intNormPcwPer == 0 && intPitchShiftPcwPer + intPitchShiftChirpPcwPer != 0){
							//if (intNormPcwPer == 0 && intPitchShiftPcwPer != 0){
							/******************** 格子シフト量チャープ(2013/7/19) ********************/

							if(y % 2 == 1){
								//								//Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#if PCW_Air_Or_SiO2
								//								Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#else
								//								Pnum[x][y].X += intPitchX; 		// intPitchXだけ+X座標に置く
								//#endif

								if (Pnum[x][y].X > 0){
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
									if (PITCH_SHIFT_MAX == PITCH){
										Pnum[x][y].X += intPitchX;
									}
									else{
										//Pnum[x][y].X += intPitchShiftX; 		// intPitchXだけ+X座標に置く
										//Pnum[x][y].X += (int) ((intPitchShiftX - intPitchX) / (double) (intPitchShiftChirpPcwPer-1)); 		// intPitchXだけ+X座標に置く
										//Pnum[x][y].X += (int) ((INT_DIV(PITCH_SHIFT_MAX2,CELL_SIZE) - intPitchX) );// (double) (PITCH_SHIFT_CHIRP_PER2-1));
									}
									/******************** 格子シフト量チャープ(2013/7/19) ********************/
								}
								else if (Pnum[x-1][y].X > 0){
									//Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}
								else if (Pnum[x-2][y].X > 0){
									//Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								}

								int poo = 0;

								if (PCW_WIDTH_CHIRP != 0){
									// PCW_WIDTH_CHIRP用微調整
									//Pnum[x][y].Y -= 0;
									poo = INT_DIV (PCW_WIDTH_CHIRP_OUT, CELL_SIZE);
									Pnum[x][y].Y -= INT_DIV(poo, (intPitchShiftPcwPerOut-1));
								}


								/******************** 格子シフト量チャープ(2013/7/19) ********************/
								if (PITCH_SHIFT_MAX != PITCH){
									//Pnum[x][y].Y = Pnum_Init[0][intPcwWid-1].Y - intPitchShiftY * (intPcwWid-1 - y) - intPCWwidthChirp;
								}
								/******************** 格子シフト量チャープ(2013/7/19) ********************/

								if (PITCH_SHIFT_MAX2 != PITCH){
									//for (y2 = intPcwWid-1; y2 >= 0; y2--){
										Pnum[x][y].X += INT_DIV(PITCH_SHIFT_MAX2,CELL_SIZE) - intPitchX;
									//}
									y_poo = 0;
									for (y2 = intPcwWid-3; y2 >= 0; y2--){
											Pnum[x][y2].Y -= INT_DIV(((PITCH_SHIFT_MAX2 * sqrt(3.0)/2) - (PITCH * sqrt(3.0)/2))/PITCH_SHIFT_CHIRP_PER2*(x - output_PitchShiftChirpPcw_Xend +1)* y_poo,CELL_SIZE) ;

										if(y2 % 2 == 1)
										y_poo++;
									}


								}


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
								//								//Pnum[x][y].X = Pnum[x-2][y].X + 3 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#if PCW_Air_Or_SiO2
								//								Pnum[x][y].X = Pnum[x-1][y].X + 2 * intPitchShiftX; 		// intPitchXだけ+X座標に置く
								//#else
								//								Pnum[x][y].X += intPitchX; 		// intPitchXだけ+X座標に置く
								//#endif
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
									//Pnum[x][y].Y -= 0;
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

					// 線欠陥から数えた偶数列円孔の端を配置する場合に使用
					//if(y % 2 == 0){
					//	if(x == 0){
					//		mcircle(Pnum[x][y].X - intPitchShiftX, Pnum[x][y].Y, z, 1);
					//	}
					//	if(x == intPcwPer - 1){
					//		mcircle(Pnum[x][y].X + intPitchShiftX, Pnum[x][y].Y, z, 1);
					//	}
					//}

				} //☆
				/****************************** LSPCW ******************************/

			}
		}
		// CORE細線導波路の位置を計算
		intWirePer2 = Pnum[intPcwPer-1][intPcwWid-1].X + intRadius;		// 出射COREスラブの開始点
		intWirePer3 = intWirePer2 + intWireLen2;						// 出射COREスラブの終了点

	}
	/****************************** フォトニック結晶 ******************************/ //739

	/****************************** 入出射細線導波路 ******************************/
	int intPcwSislabOffset;

	// 全面スラブになっているので，細線以外の部分を空気に変更
	if (PCW_SiSLAB_OFFSET != 0) {
		intPcwSislabOffset = INT_DIV(PCW_SiSLAB_OFFSET, CELL_SIZE);
	}
	else {
		intPcwSislabOffset = 0;
	}

	for (z = zmax_all - intSlabHeigPer; z < (zmax_all + 1); z++) {

		if (BOUNDARYLINE == 0) {//従来構造
								// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				if (PCW_SiSLAB_OFFSET != 0) {
					//intPcwSislabOffset;
				}
				//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset; x++){		// 配列の引数に使用するので-1
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			//出射
			for (y = 0; y < ymax_all - intWireWid_2_Out; y++) {
				/*if (PCW_SiSLAB_OFFSET != 0) {
				//intPcwSislabOffset--;		//入出射同一構造にするためのおまじない　
				}*/
				for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
					//for (x = intWirePer2 + intPcwSislabOffset+8; x < intWirePer3-8; x++){	//nondoped
					//for (x = intWirePer2 + intPcwSislabOffset-6; x < intWirePer3; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}
		}
		else if (BOUNDARYLINE == 1) {//WGテーパ円孔1列目(線形)
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - WG_chirp_in; y++) {
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x > WG_chirp_off_in_x) {
					WG_chirp_in = WG_chirp_in + WG_chirp_gradient;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < WG_chirp_off_out_x) {
					WG_chirp_out = WG_chirp_out + WG_chirp_gradient;
				}
			}

		}
		else if (BOUNDARYLINE == 2) {//WGテーパ円孔1列目(放物線)
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x)*(x - WG_chirp_off_in_x) / 17;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < (WG_chirp_off_out_x)) {
					WG_chirp_out = (-(x - WG_chirp_off_out_x)*(x - WG_chirp_off_out_x) / 17) + 16;
					//WG_chirp_out = WG_chirp_out + WG_chirp_gradient;
				}
			}

		}
		else if (BOUNDARYLINE == 3) {//WGテーパ円孔1列目(逆放物線) 2015/11/24
			double WG_chirp_in = 0;
			double WG_chirp_in_dash = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out_dash = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x)) {
					//WG_chirp_in_dash = 17 * (x - WG_chirp_off_out_y);

					WG_chirp_in_dash = 17.0*(x - WG_chirp_off_in_x);
					WG_chirp_in = -sqrt(WG_chirp_in_dash);

					//WG_chirp_in = -sqrt(100.0);
				}
			}
		}
		else if (BOUNDARYLINE == 4) {//WGテーパ円孔1列目(放物線改) 2015/11/25
			double WG_chirp_in = 0;
			double WG_chirp_in_dash = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out_dash = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x - 16)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x2 + 32)*(x - WG_chirp_off_in_x2 + 32) / 34;
				}

				/*if (x > (WG_chirp_off_in_x)) {
				WG_chirp_in = -(x - WG_chirp_off_in_x)*(x - WG_chirp_off_in_x) / 17;
				}*/ //元

			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {

				/*if (x > 0) {
				ALL_cell[x][y][z] = CORE;
				ALL_epsilonx[x][y][z] = epsilon1;
				ALL_epsilony[x][y][z] = epsilon1;
				ALL_epsilonz[x][y][z] = epsilon1;
				}*/

				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}



				if (x < (WG_chirp_off_out_x2)) {
					WG_chirp_out = (-(x - WG_chirp_off_out_x - 18)*(x - WG_chirp_off_out_x - 18) / 68) + 16;
				}


				/*if (x < (WG_chirp_off_out_x)) {
				WG_chirp_out = (-(x - WG_chirp_off_out_x)*(x - WG_chirp_off_out_x) / 17) + 16;
				//WG_chirp_out = WG_chirp_out + WG_chirp_gradient;
				}*/



			}


		}
		else if (BOUNDARYLINE == 5) {//WGテーパ円孔1列目(線形改)
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - WG_chirp_in; y++) {
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x > WG_chirp_off_in_x2) {
					WG_chirp_in = WG_chirp_in + WG_chirp_gradient2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < WG_chirp_off_out_x2) {
					WG_chirp_out = WG_chirp_out + WG_chirp_gradient2;
				}
			}

		}
		else if (BOUNDARYLINE == 6) {//WGテーパ円孔1列目(放物線改2) 2015/11/25
			double WG_chirp_in = 0;
			double WG_chirp_in_dash = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out_dash = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x - 16)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x2 + 32)*(x - WG_chirp_off_in_x2 + 32) / 68;
				}

				/*if (x > (WG_chirp_off_in_x)) {
				WG_chirp_in = -(x - WG_chirp_off_in_x)*(x - WG_chirp_off_in_x) / 17;
				}*/ //元

			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {

				/*if (x > 0) {
				ALL_cell[x][y][z] = CORE;
				ALL_epsilonx[x][y][z] = epsilon1;
				ALL_epsilony[x][y][z] = epsilon1;
				ALL_epsilonz[x][y][z] = epsilon1;
				}*/

				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}



				if (x < (WG_chirp_off_out_x2)) {
					WG_chirp_out = (-(x - WG_chirp_off_out_x - 18)*(x - WG_chirp_off_out_x - 18) / 68) + 16;
				}


				/*if (x < (WG_chirp_off_out_x)) {
				WG_chirp_out = (-(x - WG_chirp_off_out_x)*(x - WG_chirp_off_out_x) / 17) + 16;
				//WG_chirp_out = WG_chirp_out + WG_chirp_gradient;
				}*/

			}
		}
		else if (BOUNDARYLINE == 7) {//WGテーパ円孔1列目(放物線改3)
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x + 8)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x - 8)*(x - WG_chirp_off_in_x - 8) / 8;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < (WG_chirp_off_out_x)-8) {
					WG_chirp_out = (-(x - WG_chirp_off_out_x + 6)*(x - WG_chirp_off_out_x + 6) / 8) + 16;
					//WG_chirp_out = WG_chirp_out + WG_chirp_gradient;
				}
			}

		}
		else if (BOUNDARYLINE == 8) {//國分先生の提案構造
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x + 4)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x - 4)*(x - WG_chirp_off_in_x - 4) / 32;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < (WG_chirp_off_out_x)-4) {
					WG_chirp_out = (-(x - WG_chirp_off_out_x + 4)*(x - WG_chirp_off_out_x + 4) / 32) + 16;
					//WG_chirp_out = WG_chirp_out + WG_chirp_gradient;
				}
			}

		}
		else if (BOUNDARYLINE == 9) { //WG緩やかテーパ1(出射のみ)
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < WG_chirp_off_out_x2 + 50) {
					WG_chirp_out = WG_chirp_out + WG_chirp_gradient3;
				}
			}

		}
		else if (BOUNDARYLINE == 10) { //WGさらに緩やかテーパ1(出射のみ)
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x < WG_chirp_off_out_x2 + 54) {
					WG_chirp_out = WG_chirp_out + WG_chirp_gradient4;
				}
			}

		}
		else if (BOUNDARYLINE == 11) { //放物線兼さらにWG緩やかテーパ(出射のみ)2015/12/15
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x <= WG_chirp_off_out_x2 - 19) {
					WG_chirp_out = 6;
				}

				if (WG_chirp_off_out_x2 - 19 <= x) { //876
					if (x <= WG_chirp_off_out_x2 + 55) {

						if (x == WG_chirp_off_out_x2 - 19) {
							//WG_chirp_out = 0;
						}

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient4;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 12) { //PCW近接オフセット深く兼さらにWG緩やかテーパ(出射のみ)2015/12/15
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x <= WG_chirp_off_out_x2 - 18) {
					WG_chirp_out = 3;
				}

				if (x == WG_chirp_off_out_x2 - 17 - 12 - 3) {
					WG_chirp_out = 2;
				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x2 + 55) {

						if (x == WG_chirp_off_out_x2 - 17 - 12) {
							WG_chirp_out = 0;
						}

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient4;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 13) { //12+出射放物線2015/12/15
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x + 8)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x - 8)*(x - WG_chirp_off_in_x - 8) / 8;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x <= WG_chirp_off_out_x2 - 18) {
					WG_chirp_out = 3;
				}

				if (x == WG_chirp_off_out_x2 - 17 - 12 - 3) {
					WG_chirp_out = 2;
				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x2 + 55) {

						if (x == WG_chirp_off_out_x2 - 17 - 12) {
							WG_chirp_out = 0;
						}

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient4;
					}
				}
			}//出射終


		}
		else if (BOUNDARYLINE == 14) { //PCW近接オフセット深く兼さらにWG緩やかテーパ(出射のみ) + 出射モニター後ろ2015/12/15
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9 + 4 + 1; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				/*if (x <= WG_chirp_off_out_x2 - 18) {
				WG_chirp_out = 3;
				}*/

				/*if (x == WG_chirp_off_out_x2 - 17 - 12 - 3) {
				WG_chirp_out = 2;
				}*/

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x9) {

						/*if (x == WG_chirp_off_out_x9 - 17 - 12) {
						WG_chirp_out = 0;
						}*/

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient5;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 141) { //PCW近接オフセット深く兼さらにWG緩やかテーパ(出射のみ) + 出射モニター後ろ2015/12/15
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y + 9 + 7; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}


			}//出射終

		}
		else if (BOUNDARYLINE == 15) { //テーパ変化2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y15 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x9) {

						/*if (x == WG_chirp_off_out_x9 - 17 - 12) {
						WG_chirp_out = 0;
						}*/

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient6;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 16) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y16 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x9) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient7;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 17) { //テーパ変化 2015/12/24
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y17 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x9) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient8;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 18) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y18 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 - 100 <= x) { //876
					if (x <= WG_chirp_off_out_x18) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient9;
					}
				}
			}//出射終

		}
		else if (BOUNDARYLINE == 19) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y19 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x19) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient10;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 20) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient11;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 21) { //入射テーパ 2016/1/4
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 + WG_chirp_in; y++) {

					//for (y = 0; y < ymax_all - intWireWid_2; y++) {
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in; x++) {		// 配列の引数に使用するので-1
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped


					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
				if (x >(WG_chirp_off_in_x21)) {
					WG_chirp_in = -(x - WG_chirp_off_in_x21)*(x - WG_chirp_off_in_x21) / 8;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y21 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x21) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient11;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 22) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient12;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 23) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient13;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 24) { //テーパ変化 2015/12/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient14;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 25) { //テーパ変化 1/7
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient15;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 26) { //テーパ変化 1/7
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient16;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 27) { //テーパ変化 1/7
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4 + 2; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {


						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 28) { //テーパなし2/5
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2 - 10; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][ y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 + 9 + 4 + 2; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 29) { //テーパなし2/5
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2 - 10 - 2; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 30) { //テーパなし2/5
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2 - 10 - 2 - 5; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 31) { //テーパなし2/8
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2 - 10 - 2 - 5-10; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5-10; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 32) { //テーパなし+ 出射テーパ(21)2/16
			double WG_chirp_out = 0;

			// 入射
			for (y = 0; y < ymax_all - intWireWid_2 - 10 - 2 - 5 - 10; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																						//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient11;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 33) { //テーパなし+ 入射テーパ3/3
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 10 - 2 - 5 - 10 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (x >(WG_chirp_off_in_x33 + 8)) {
				 WG_chirp_in = -(x - WG_chirp_off_in_x33 - 8)*(x - WG_chirp_off_in_x33 - 8) / 8;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}else if (BOUNDARYLINE == 34) { //SiWGテーパ 傾き0.08
	double WG_chirp_in = 0;
	double WG_chirp_out = 0;

	// 入射
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in; y++) {
			// 配列の引数に使用するので-1
			//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
			//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
		if (x > WG_chirp_off_in_x34-3) {
			WG_chirp_in = WG_chirp_in + 0.08;
		}
	}

	// 出射
	for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
		for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - 20 /*+WG_chirp_off_out_y20*/; y++) {//16/1/26

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;

		}

		if (WG_chirp_off_out_x2 + 290 - 30 - 80 - 15-85 <= x) { //876
			if (WG_chirp_off_out_x2 + 270-70-10 > x) { //876
				WG_chirp_out = WG_chirp_out - 0.08;
			}
		}
	}//出射終
}else if (BOUNDARYLINE == 35) { //SiWGテーパ 傾き0.06
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in+2; y++) {
					// 配列の引数に使用するので-1
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
				if (x > WG_chirp_off_in_x34+1) {
					WG_chirp_in = WG_chirp_in + 0.06;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - 20 /*+WG_chirp_off_out_y20*/; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 + 290 - 30 - 80 - 15 - 85 <= x) { //876
					if (WG_chirp_off_out_x2 + 270 - 70-10 > x) { //876
						WG_chirp_out = WG_chirp_out - 0.06;
					}
				}
			}//出射終
		}else if (BOUNDARYLINE == 36) { //SiWGテーパ 傾き0.04
			double WG_chirp_in = 0;
			double WG_chirp_out = 0;

			// 入射
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in+3; y++) {
					// 配列の引数に使用するので-1
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
					//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
				if (x > WG_chirp_off_in_x34 +1-17) {
					WG_chirp_in = WG_chirp_in + 0.04;
				}
			}

			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - 20 /*+WG_chirp_off_out_y20*/; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}

				if (WG_chirp_off_out_x2 + 290 - 30 - 80 - 15 - 85 <= x) { //876
					if (WG_chirp_off_out_x2 + 270 - 70 - 0 > x) { //876
						WG_chirp_out = WG_chirp_out - 0.04;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 37) { //SiWGテーパ
			double WG_chirp_in = 0;
			double WG_chirp_in2 = 0;
			double WG_chirp_out = 0;

			// 入射
			//土台部
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in ; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				/*if (x > WG_chirp_off_in_x34 -60) {
					WG_chirp_in = WG_chirp_in + 0.1;
				}*/
			}

			//可変部
		for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 4; x++) { //PCW側の端の座標

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				if (x > 300) {
					if (y > 97) {
						WG_chirp_in2 = WG_chirp_in2 + 1.75;
					}
							}
					}



			// 入射の円孔を消すプログラム
			for (x = 0; x <431; x++) {
				for (y = 0; y < ymax_all; y++) {

					if (y>155) {
						ALL_cell[x][y][z] = CORE;
						ALL_epsilonx[x][y][z] = epsilon1;
						ALL_epsilony[x][y][z] = epsilon1;
						ALL_epsilonz[x][y][z] = epsilon1;
					}
				}
			}


			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}
		else if (BOUNDARYLINE == 38) { //SiWGテーパ 幅テーパ168nm

			double WG_chirp_in = 0;
			double WG_chirp_in2 = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out2 = 0;

			// 入射
			//土台部
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				/*if (x > WG_chirp_off_in_x34 -60) {
				WG_chirp_in = WG_chirp_in + 0.1;
				}*/
			}

			//可変部
			for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 4; x++) { //PCW側の端の座標

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				if (x > 300) {
					if (y > 97) {
						WG_chirp_in2 = WG_chirp_in2 + 1.75;
					}
				}
			}



			// 入射の円孔を消すプログラム
			for (x = 0; x <431; x++) {
				for (y = 0; y < ymax_all; y++) {

					if (y>147) {
						ALL_cell[x][y][z] = CORE;
						ALL_epsilonx[x][y][z] = epsilon1;
						ALL_epsilony[x][y][z] = epsilon1;
						ALL_epsilonz[x][y][z] = epsilon1;
					}
				}
			}


			// 出射
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {//16/1/26

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;

				}
				if (WG_chirp_off_out_x2 - 17 <= x) { //876
					if (x <= WG_chirp_off_out_x20) {

						WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
					}
				}
			}//出射終
		}else if (BOUNDARYLINE == 39) { //SiWGテーパ入出射
			double WG_chirp_in = 0;
			double WG_chirp_in2 = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out2 = 0;

			// 入射
			//土台部
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				/*if (x > WG_chirp_off_in_x34 -60) {
				WG_chirp_in = WG_chirp_in + 0.1;
				}*/
			}

			//可変部
			for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 4; x++) { //PCW側の端の座標

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				if (x > 300) {
					if (y > 97) {
						WG_chirp_in2 = WG_chirp_in2 + 1.75;
					}
				}
			}



			// 入射の円孔を消すプログラム
			for (x = 0; x <431; x++) {
				for (y = 0; y < ymax_all; y++) {

					if (y>155) {
						ALL_cell[x][y][z] = CORE;
						ALL_epsilonx[x][y][z] = epsilon1;
						ALL_epsilony[x][y][z] = epsilon1;
						ALL_epsilonz[x][y][z] = epsilon1;
					}
				}
			}


			// 出射 土台
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
					}
				}//出射終

			// 出射可変部
			for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
			  for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2; x < intWirePer3; x++) {


					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				if (x > 818) {
					if (y > 94) {
						WG_chirp_out2 = WG_chirp_out2 - 1.75;
					}
				}
			}
			//出射終


			// 出射の円孔を消すプログラム
			for (x = 823; x <xmax_all; x++) {
				for (y = 0; y < ymax_all; y++) {

					if (y>155) {
						ALL_cell[x][y][z] = CORE;
						ALL_epsilonx[x][y][z] = epsilon1;
						ALL_epsilony[x][y][z] = epsilon1;
						ALL_epsilonz[x][y][z] = epsilon1;
					}
				}
			}

		}else if (BOUNDARYLINE == 40) { //SiWGテーパ入出射
			double WG_chirp_in = 0;
			double WG_chirp_in2 = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out2 = 0;

			// 入射
			//土台部
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				/*if (x > WG_chirp_off_in_x34 -60) {
				WG_chirp_in = WG_chirp_in + 0.1;
				}*/
			}

			//可変部
			for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 4; x++) { //PCW側の端の座標

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				if (x > 300) {
					if (y > 97) {
						WG_chirp_in2 = WG_chirp_in2 + 1.;
					}
				}
			}



			// 入射の円孔を消すプログラム
			for (x = 0; x <390; x++) {
				for (y = 0; y < ymax_all; y++) {

					if (y>155) {
						ALL_cell[x][y][z] = CORE;
						ALL_epsilonx[x][y][z] = epsilon1;
						ALL_epsilony[x][y][z] = epsilon1;
						ALL_epsilonz[x][y][z] = epsilon1;
					}
				}
			}


			// 出射 土台
			for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
				for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}//出射終

			 // 出射可変部
			for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
				for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2; x < intWirePer3; x++) {


					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				if (x > 818) {
					if (y > 94) {
						WG_chirp_out2 = WG_chirp_out2 - 1.3;
					}
				}
			}
			//出射終


			// 出射の円孔を消すプログラム
			for (x = 840; x <xmax_all; x++) {
				for (y = 0; y < ymax_all; y++) {

					if (y>155) {
						ALL_cell[x][y][z] = CORE;
						ALL_epsilonx[x][y][z] = epsilon1;
						ALL_epsilony[x][y][z] = epsilon1;
						ALL_epsilonz[x][y][z] = epsilon1;
					}
				}
			}

		}
		else if (BOUNDARYLINE == 41) { //SiWGテーパ入出射 円孔調整有 外側接触
			double WG_chirp_in = 0;
			double WG_chirp_in2 = 0;
			double WG_chirp_out = 0;
			double WG_chirp_out2 = 0;

			// 入射
			//土台部
			for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-100; x++) {
				for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in+5; y++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

				//可変部
				for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
					for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1+5; x++) { //PCW側の端の座標

						ALL_cell[x][y][z] = CLAD;
						ALL_epsilonx[x][y][z] = epsilon2;
						ALL_epsilony[x][y][z] = epsilon2;
						ALL_epsilonz[x][y][z] = epsilon2;
					}

					//if (x > 220) {
					if (y > 79) //テーパ開始のy座標
						WG_chirp_in2 = WG_chirp_in2 + 1.33;
						//}
					}

			//可変部
			for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
				for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1; x++) { //PCW側の端の座標

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}

			 // 出射可変部
			for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
				for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6+6-20+3; x < intWirePer3-20; x++) {


					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}

				//if (x > 818) {
					if (y > 94-15) {
						WG_chirp_out2 = WG_chirp_out2 - 1.33;
					}
				//}
			}
			//可変部
			for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
				for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3+5; x < intWirePer3 - 20; x++) {

					ALL_cell[x][y][z] = CLAD;
					ALL_epsilonx[x][y][z] = epsilon2;
					ALL_epsilony[x][y][z] = epsilon2;
					ALL_epsilonz[x][y][z] = epsilon2;
				}
			}



			//// 出射の円孔を消すプログラム //これは使う必要はない
			//for (x = 823; x <xmax_all; x++) {
			//	for (y = 0; y < ymax_all; y++) {

			//		if (y>155) {
			//			ALL_cell[x][y][z] = CORE;
			//			ALL_epsilonx[x][y][z] = epsilon1;
			//			ALL_epsilony[x][y][z] = epsilon1;
			//			ALL_epsilonz[x][y][z] = epsilon1;
			//		}
			//	}
			//}

		}
else if (BOUNDARYLINE == 42) { //SiWGテーパ入出射 円孔調整有 中心接触
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5-5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.33;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22-5; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1-5; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}



	// 入射の円孔を消すプログラム
	for (x = 0; x <431; x++) {
		for (y = 0; y < ymax_all; y++) {

			if (y>155) {
				ALL_cell[x][y][z] = CORE;
				ALL_epsilonx[x][y][z] = epsilon1;
				ALL_epsilony[x][y][z] = epsilon1;
				ALL_epsilonz[x][y][z] = epsilon1;
			}
		}
	}


	// 出射 土台
	for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
		for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10-5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}//出射終

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3; x < intWirePer3 - 20; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.33;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22-5; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}



	// 出射の円孔を消すプログラム
	for (x = 823; x <xmax_all; x++) {
		for (y = 0; y < ymax_all; y++) {

			if (y>155) {
				ALL_cell[x][y][z] = CORE;
				ALL_epsilonx[x][y][z] = epsilon1;
				ALL_epsilony[x][y][z] = epsilon1;
				ALL_epsilonz[x][y][z] = epsilon1;
			}
		}
	}

}
else if (BOUNDARYLINE == 43) { //SiWGテーパ入出射 円孔調整有 外側接触
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 =0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5-5-5; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.33;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10+2; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}



	//// 入射の円孔を消すプログラム //これは使う必要はない
	//for (x = 0; x <431; x++) {
	//	for (y = 0; y < ymax_all; y++) {

	//		if (y>155) {
	//			ALL_cell[x][y][z] = CORE;
	//			ALL_epsilonx[x][y][z] = epsilon1;
	//			ALL_epsilony[x][y][z] = epsilon1;
	//			ALL_epsilonz[x][y][z] = epsilon1;
	//		}
	//	}
	//}


	//// 出射 土台
	//for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
	//	for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {

	//		ALL_cell[x][y][z] = CLAD;
	//		ALL_epsilonx[x][y][z] = epsilon2;
	//		ALL_epsilony[x][y][z] = epsilon2;
	//		ALL_epsilonz[x][y][z] = epsilon2;
	//	}
	//}//出射終

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3+10; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.33;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5+10-3; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}


}
else if (BOUNDARYLINE == 44) { //SiWGテーパ入出射 円孔調整有 外側接触 2周期
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5+10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 95) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.55;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2+9; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 10-7-3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 95) {
			WG_chirp_out2 = WG_chirp_out2 - 1.55;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 -9-3; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}
}else if (BOUNDARYLINE == 45) { //SiWGテーパ入出射 円孔調整有 外側接触 1.5周期
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 95+16) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 2.09;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2 + 9+1; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 10 - 7 - 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 95+16) {
			WG_chirp_out2 = WG_chirp_out2 - 2.09;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 9 - 3-1; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}
}
else if (BOUNDARYLINE == 46) { //SiWGテーパ入出射 円孔調整有 外側接触 1周期
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 95 + 16) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 2.09;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2 + 9 + 1; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 10 - 7 - 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 95 + 16) {
			WG_chirp_out2 = WG_chirp_out2 - 2.09;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 9 - 3 - 1; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}
}
else if (BOUNDARYLINE == 47) { //SiWGテーパ入出射 円孔調整有 外側接触(2.5周期) x=8a
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5+10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.45;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2+10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.45;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 3-10; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}


}
else if (BOUNDARYLINE == 48) { //SiWGテーパ入出射 円孔調整有 外側接触(2.5周期) x=9a
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.75;
		//}
	}

	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2 + 10-3; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.75;
		}
		//}
	}
	//可変部
	for (y = 155; y < ymax_all - intWireWid_2 - 22; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 3 - 10-3+6; x < intWirePer3 - 20-7; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}


}
else if (BOUNDARYLINE == 49) { //テーパなし11/4
	double WG_chirp_in = 0;
	double WG_chirp_out = 0;

	// 入射
	for (y = 0; y < ymax_all - intWireWid_2 - 10 - 2 - 5 - 10; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++) {		// 配列の引数に使用するので-1
																				//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1-8; x++){		//nondoped
																				//for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1+6; x++){		//doped
			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 22-5; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in - 10 + 2 + 10 - 3+579; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	// 出射
	for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
		for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {//16/1/26

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;

		}

		if (WG_chirp_off_out_x2 - 17 <= x) { //876
			if (x <= WG_chirp_off_out_x20) {

				WG_chirp_out = WG_chirp_out + WG_chirp_gradient17;
			}
		}
	}//出射終
}
else if (BOUNDARYLINE == 50) { //SiWGテーパ入出射 円孔調整有 外側接触 幅チャ有 N =7
//PCW1列目はwtaper 21 nmに相当
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27+5; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5+10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.25;
		//}
	}

	//可変部
	for (y = 153; y < ymax_all - intWireWid_2 - 22+6; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2+10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}



	//// 入射の円孔を消すプログラム //これは使う必要はない
	//for (x = 0; x <431; x++) {
	//	for (y = 0; y < ymax_all; y++) {

	//		if (y>155) {
	//			ALL_cell[x][y][z] = CORE;
	//			ALL_epsilonx[x][y][z] = epsilon1;
	//			ALL_epsilony[x][y][z] = epsilon1;
	//			ALL_epsilonz[x][y][z] = epsilon1;
	//		}
	//	}
	//}


	//// 出射 土台
	//for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
	//	for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {

	//		ALL_cell[x][y][z] = CLAD;
	//		ALL_epsilonx[x][y][z] = epsilon2;
	//		ALL_epsilony[x][y][z] = epsilon2;
	//		ALL_epsilonz[x][y][z] = epsilon2;
	//	}
	//}//出射終

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10+5; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.25;
		}
		//}
	}
	//可変部
	for (y = 153; y < ymax_all - intWireWid_2 - 22+6; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 3-10; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}


}
else if (BOUNDARYLINE == 51) { //SiWGテーパ入出射 円孔調整有 外側接触 幅チャ有 N = 9 テーパ周期を6箇所全て変更 (LSPCW_PERは-4して12に変更)
//PCW1列目はwtaper 63 nmに相当
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27 + 5; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.25;
		//}
	}

	//可変部
	for (y = 153; y < ymax_all - intWireWid_2 - 22 + 6-1; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}
	//// 入射の円孔を消すプログラム //これは使う必要はない
	//for (x = 0; x <431; x++) {
	//	for (y = 0; y < ymax_all; y++) {

	//		if (y>155) {
	//			ALL_cell[x][y][z] = CORE;
	//			ALL_epsilonx[x][y][z] = epsilon1;
	//			ALL_epsilony[x][y][z] = epsilon1;
	//			ALL_epsilonz[x][y][z] = epsilon1;
	//		}
	//	}
	//}


	//// 出射 土台
	//for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
	//	for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {

	//		ALL_cell[x][y][z] = CLAD;
	//		ALL_epsilonx[x][y][z] = epsilon2;
	//		ALL_epsilony[x][y][z] = epsilon2;
	//		ALL_epsilonz[x][y][z] = epsilon2;
	//	}
	//}//出射終

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10 + 5; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.25;
		}
		//}
	}
	//可変部
	for (y = 153; y < ymax_all - intWireWid_2 - 22 + 6-1; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 3 - 10; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}


}
else if (BOUNDARYLINE == 52) { //SiWGテーパ入出射 円孔調整有 外側接触 幅チャ有 N = 9 テーパ周期を6箇所全て変更 (LSPCW_PERは-4して12に変更)
							   //PCW1列目はwtaper 63 nmに相当 幅チャは途中から
	double WG_chirp_in = 0;
	double WG_chirp_in2 = 0;
	double WG_chirp_out = 0;
	double WG_chirp_out2 = 0;

	// 入射
	//土台部
	for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset - 1 - 100; x++) {
		for (y = 0; y < ymax_all - intWireWid_2 - 27 + WG_chirp_in + 5; y++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}

	//可変部
	for (y = 0; y < ymax_all - intWireWid_2 - 27 + 5; y++) {
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 1 + 5 - 5 - 5 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 220) {
		if (y > 79) //テーパ開始のy座標
			WG_chirp_in2 = WG_chirp_in2 + 1.25;
		//}
	}

	//可変部
	for (y = 153; y < ymax_all - intWireWid_2 - 22 + 6 - 1; y++) {//長法領域
		for (x = 0; x < intWireLen1 - 1 - intPcwSislabOffset + WG_chirp_in2 - 10 + 2 + 10; x++) { //PCW側の端の座標

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}
	//// 入射の円孔を消すプログラム //これは使う必要はない
	//for (x = 0; x <431; x++) {
	//	for (y = 0; y < ymax_all; y++) {

	//		if (y>155) {
	//			ALL_cell[x][y][z] = CORE;
	//			ALL_epsilonx[x][y][z] = epsilon1;
	//			ALL_epsilony[x][y][z] = epsilon1;
	//			ALL_epsilonz[x][y][z] = epsilon1;
	//		}
	//	}
	//}


	//// 出射 土台
	//for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer3; x++) {
	//	for (y = 0; y < ymax_all - intWireWid_2_Out + WG_chirp_out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10; y++) {

	//		ALL_cell[x][y][z] = CLAD;
	//		ALL_epsilonx[x][y][z] = epsilon2;
	//		ALL_epsilony[x][y][z] = epsilon2;
	//		ALL_epsilonz[x][y][z] = epsilon2;
	//	}
	//}//出射終

	// 出射可変部
	for (y = 0; y < ymax_all - intWireWid_2_Out - WG_chirp_off_out_y20 - 2 + 9 + 4 + 2 - 5 - 10 + 5; y++) {
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3; x < intWirePer3; x++) {


			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}

		//if (x > 818) {
		if (y > 94 - 15) {
			WG_chirp_out2 = WG_chirp_out2 - 1.25;
		}
		//}
	}
	//可変部
	for (y = 153; y < ymax_all - intWireWid_2 - 22 + 6 - 1; y++) {//長法領域
		for (x = intWirePer2 + intPcwSislabOffset + WG_chirp_out2 + 6 + 6 - 20 + 3 + 5 + 10 - 3 - 10; x < intWirePer3 - 20; x++) {

			ALL_cell[x][y][z] = CLAD;
			ALL_epsilonx[x][y][z] = epsilon2;
			ALL_epsilony[x][y][z] = epsilon2;
			ALL_epsilonz[x][y][z] = epsilon2;
		}
	}


}




		//田村
		//入射
		/*for (y = 0; y < ymax_all - intWireWid_2; y++){
		for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -8; x < intWireLen1 - 1 - intPcwSislabOffset - 7; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-1; y++){
		for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -6; x < intWireLen1 - 1 - intPcwSislabOffset - 5; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-2; y++){
		for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -4; x < intWireLen1 - 1 - intPcwSislabOffset - 3; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-4; y++){
		for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -2; x < intWireLen1 - 1 - intPcwSislabOffset - 2; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
		for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -1; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
		for (x = intWireLen1 - 1 - intPcwSislabOffset - 1 -1; x < intWireLen1 - 1 - intPcwSislabOffset - 1; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		//出射
		for (y = 0; y < ymax_all - intWireWid_2; y++){
		for (x = intWirePer2 + intPcwSislabOffset+6; x < intWirePer2 + intPcwSislabOffset+8; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-1; y++){
		for (x = intWirePer2 + intPcwSislabOffset+4; x < intWirePer2 + intPcwSislabOffset+6; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-2; y++){
		for (x = intWirePer2 + intPcwSislabOffset+2; x < intWirePer2 + intPcwSislabOffset+4; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-4; y++){
		for (x = intWirePer2 + intPcwSislabOffset+1; x < intWirePer2 + intPcwSislabOffset+2; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
		for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer2 + intPcwSislabOffset+1; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}
		for (y = 0; y < ymax_all - intWireWid_2-6; y++){
		for (x = intWirePer2 + intPcwSislabOffset; x < intWirePer2 + intPcwSislabOffset+1; x++){		// 配列の引数に使用するので-1
		ALL_cell[x][y][z] = CLAD;
		ALL_epsilonx[x][y][z] = epsilon2;
		ALL_epsilony[x][y][z] = epsilon2;
		ALL_epsilonz[x][y][z] = epsilon2;
		}
		}*/
		//田村ここまで
	}
	/****************************** 入出射細線導波路 ******************************/ //1915



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
	/****************************** 各ノードにモデルを分割 ******************************/

	/****************************** 共通パラメータの設定 ******************************/

	// 励振点，観測面の設定 (XMAXは "のりしろ" 部分を含めていることに注意)
	intExctPortNum = intExctLen / (XMAX - 1);
	intObseInPortNum = intObseLen1 / (XMAX - 1);
	intObseOutPortNum = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) / (XMAX - 1);
	if (NODE % 2 != 0) {
		intObseCenPortNum = XMAX_ALL / 2 / (XMAX - 1);	// 奇数個並列計算のとき
	}
	else {
		intObseCenPortNum = XMAX_ALL / 2 / (XMAX - 1) - 1;	// 偶数個並列計算のとき
	}

	intExctLenPart = intExctLen % (XMAX - 1) - 1;		// 配列の引数に使用するので-1 //☆入射モニタ
	intObseLenPart1 = intObseLen1 % (XMAX - 1) - INT_DIV(intObseInter, 2) - 1;
	intObseLenPart2 = intObseLenPart1 + intObseInter;
	intObseLenPart3 = intObseLen1 % (XMAX - 1);
	intObseLenPart4 = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1) - INT_DIV(intObseInter, 2);
	//☆☆出射モニタが移動．これだ．	+ 119(セルサイズ)　+ で出射端側に移動
	//× ここをいじるとずれるのでやらない

	intObseLenPart5 = intObseLenPart4 + intObseInter;
	intObseLenPart6 = (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1);//☆中点とは？→おそらく真中にある緑色の線のこと

																					 //intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1) - INT_DIV(intObseInter, 2);
																					 //intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1);

																					 //intObseLenPart10  =  (intWirePer2 + INT_DIV(OBSE_WIRE_LEN, CELL_SIZE)) % (XMAX - 1)  - 20;//☆☆PCW内のモニタ 16/1/5


	if (NODE % 2 != 0) {
		intObseLenPart7 = (XMAX_ALL / 2) % (XMAX - 1) - 10;		// 奇数個並列計算のとき
	}
	else {
		intObseLenPart7 = (XMAX - 1) - 10;		// 偶数個並列計算のとき
	}
	//intObseLenPart8 = intObseLenPart7 + intObseInter;
	//intObseLenPart9 = (XMAX_ALL / 2) % (XMAX - 1);

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
		for(z = 0; z < zmax+1; z++){
			epsilon_yz[y][z] = epsilony[intObseLenPart1][y][z];
		}
	}
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilony, "%e\t", epsilon_yz[y][z]);
		}
		fprintf(fpepsilony, "\n");
	}

	//ZX平面 (Y:境界面)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx[x][z] = epsilonz[x][ymax][z];
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilonz, "%e\t", epsilon_zx[x][z]);
		}
		fprintf(fpepsilonz, "\n");
	}
	//ZX平面 (Y:中心)
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			epsilon_zx2[x][z] = epsilonz[x][ymax/2][z];
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			fprintf(fpepsilonz2, "%e\t", epsilon_zx2[x][z]);
		}
		fprintf(fpepsilonz2, "\n");
	}

	////YZ平面
	//for(y = 0; y < ymax+1; y++){
	//	for(z = 0; z <= zmax; z++){
	//		epsilon_yz[y][z] = epsilony[intObseLenPart1][y][z];
	//	}
	//}
	//for(y = 0; y < ymax+1; y++){
	//	for(z = 0; z <= zmax; z++){
	//		fprintf(fpepsilony, "%e\t", epsilon_yz[y][z]);
	//	}
	//	fprintf(fpepsilony, "\n");
	//}

	////ZX平面 (Y:境界面)
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		epsilon_zx[x][z] = epsilonz[x][ymax][z];
	//	}
	//}
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		fprintf(fpepsilonz, "%e\t", epsilon_zx[x][z]);
	//	}
	//	fprintf(fpepsilonz, "\n");
	//}
	////ZX平面 (Y:中心)
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		epsilon_zx2[x][z] = epsilonz[x][ymax/2][z];
	//	}
	//}
	//for(x = 0; x < xmax; x++){
	//	for(z = 0; z < zmax+1; z++){
	//		fprintf(fpepsilonz2, "%e\t", epsilon_zx2[x][z]);
	//	}
	//	fprintf(fpepsilonz2, "\n");
	//}


	// ファイルポインタを閉じる
	if (irank == IRANK_MIN){
		fclose(fpallepsilonx);
	}
	fclose(fpepsilonx);
	fclose(fpepsilony);
	fclose(fpepsilonz);


	//MPI_Barrier(MPI_COMM_WORLD); 			// 一度同期をとる

	//// メモリの開放
	//for(x = 0; x < XMAX+1; x++)	free(epsilon_xy[x]);
	//for(y = 0; y < YMAX+1; y++)	free(epsilon_yz[y]);
	//for(x = 0; x < XMAX+1; x++)	free(epsilon_zx[x]);

	//free(epsilon_xy);
	//free(epsilon_yz);
	//free(epsilon_zx);
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
		for(z = 0; z < zmax; z++){
			Hx[x][ymax][z] = Hx[x][ymax-1][z];		// 偶関数
		}
	}
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			Hz[x][ymax][z] = Hz[x][ymax-1][z];		// 偶関数
		}
	}
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			Hy[x][y][zmax] = -Hy[x][y][zmax-1];		// 奇関数
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Hx[x][y][zmax] = -Hx[x][y][zmax-1];		// 奇関数
		}
	}

	/****************************** 磁界の対称境界条件(4回対称) ******************************/
}


// モデルへの励振点，観測点の記録
void observation_func(){

	//int y, z;

	if(irank == intObseInPortNum){ //入射

		for(int x = intObseLenPart1; x < intObseLenPart2; x++){
			/****************************** 観測面の修正(2013/8/8) ******************************/
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
					//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
					/****************************** 観測面の修正(2013/8/8) ******************************/
					if((y == YMAX-1) && (z == (intSlabCen-1))){		//活性層断面中央点の時を記録
						cell[x][y][z] = 4; 					//中央点確認用
					}
					cell[x][y][z] += OBSERVATION; 		//面確認用
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
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
					//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
					/****************************** 観測面の修正(2013/8/8) ******************************/

					if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
						cell[x][y][z] = 4; 					// 中央点確認用
					}
					cell[x][y][z] += OBSERVATION; 			// 面確認用
				}
			}
		}
	}

	if(irank == intObseCenPortNum){ //出射 NODE 2
		//for(int x = intObseLenPart7; x < intObseLenPart8; x++){
		int x = intObseLenPart7;
		/****************************** 観測面の修正(2013/8/8) ******************************/
		for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
			for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
				//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
				//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
				/****************************** 観測面の修正(2013/8/8) ******************************/

				if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
					cell[x][y][z] = 4; 					// 中央点確認用
				}
				cell[x][y][z] += OBSERVATION; 			// 面確認用
			}
		}

	}
}


/*if (irank == intObseCenPortNum) { //出射 NODE 2　//16/1/5
int x = intObseLenPart10;*/
/****************************** 観測面の修正(2013/8/8) ******************************/
/*for (int y = ymax - intObseWid; y < ymax; y++) { // 矩形導波路断面Y領域判断．
for (int z = zmax - intObseHeig; z < zmax; z++) {		//矩形導波路断面Z領域判断．
//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
/****************************** 観測面の修正(2013/8/8) ******************************/
/*
if ((y == YMAX - 1) && (z == (intSlabCen - 1))) {		// 活性層断面中央点の時を記録
cell[x][y][z] = 4; 					// 中央点確認用
}
cell[x][y][z] += OBSERVATION; 			// 面確認用
}
}

}

*/


void calc_efield(){

	double dex, dey, dez;
	double cnstEx, cnstEy, cnstEz;

	// Ex
	for(x = 0; x < xmax; x++){
		for(y = 1; y < ymax+1; y++){		// Exはy軸に対して奇関数
			for(z = 1; z < zmax+1; z++){
				cnstEx = dt / epsilonx[x][y][z];
				dex = ( (Hz[x][y][z] - Hz[x][y-1][z]) / dy) - ( (Hy[x][y][z] - Hy[x][y][z-1]) / dz);
				Ex[x][y][z] = Ex[x][y][z] + cnstEx * dex;
			}
		}
	}

	// Ey
	for(x = 1; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 1; z < zmax+1; z++){
				cnstEy = dt / epsilony[x][y][z];
				dey = ( (Hx[x][y][z] - Hx[x][y][z-1]) / dz)-( (Hz[x][y][z] - Hz[x-1][y][z]) / dx);
				Ey[x][y][z] = Ey[x][y][z] + cnstEy * dey;
			}
		}
	}

	// Ez
	for(x = 1; x < xmax; x++){
		for(y = 1; y < ymax+1; y++){		// Ezはy軸に対して奇関数
			for(z = 0; z < zmax; z++){
				cnstEz = dt / epsilonz[x][y][z];
				dez = ( (Hy[x][y][z] - Hy[x-1][y][z]) / dx) - ( (Hx[x][y][z] - Hx[x][y-1][z]) / dy);
				Ez[x][y][z] = Ez[x][y][z] + cnstEz * dez;
			}
		}
	}

		/*if (irank == intObseCenPortNum) { //出射 NODE 2　//16/1/5
		int x = intObseLenPart10;*/
		/****************************** 観測面の修正(2013/8/8) ******************************/
		/*for (int y = ymax - intObseWid; y < ymax; y++) { // 矩形導波路断面Y領域判断．
			for (int z = zmax - intObseHeig; z < zmax; z++) {		//矩形導波路断面Z領域判断．
																	//for(int y = 0; y <= YMAX-1; y++){ //矩形導波路断面Y領域判断．
																	//	for(int z = (air_hc+intCladHeight1); z <= (air_hc+intCladHeight1+intSlabHeigPer); z++){		//矩形導波路断面Z領域判断．-1は配列が0開始なため
																	/****************************** 観測面の修正(2013/8/8) ******************************/
	/*
				if ((y == YMAX - 1) && (z == (intSlabCen - 1))) {		// 活性層断面中央点の時を記録
					cell[x][y][z] = 4; 					// 中央点確認用
				}
				cell[x][y][z] += OBSERVATION; 			// 面確認用
			}
		}

	}

	*/


	/****************************** 電界の対称境界条件 ******************************/

	// 境界面で反対称となる電界成分の境界面上の値を0としている
	for(x = 0; x < xmax; x++){
		for(z = 0; z < zmax+1; z++){
			Ex[x][ymax][z] = 0.0;		// 奇関数
		}
	}
	for(x = 0; x < xmax+1; x++){
		for(z = 0; z < zmax; z++){
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
			for(z = 0; z < zmax; z++){
				dhx = ( (Ey[x][y][z+1] - Ey[x][y][z]) / dz) - ( (Ez[x][y+1][z] - Ez[x][y][z]) / dy);
				Hx[x][y][z] = Hx[x][y][z] + cnstHxyz * dhx;
			}
		}
	}

	// Hy
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax+1; y++){
			for(z = 0; z < zmax; z++){
				dhy = ( (Ez[x+1][y][z] - Ez[x][y][z]) / dx) - ( (Ex[x][y][z+1] - Ex[x][y][z]) / dz);
				Hy[x][y][z] = Hy[x][y][z] + cnstHxyz * dhy;
			}
		}
	}

	// Hz
	for(x = 0; x < xmax; x++){
		for(y = 0; y < ymax; y++){
			for(z = 0; z < zmax+1; z++){
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
	for(z = 0; z < zmax+1; z++){
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
	for(z = 0; z < zmax; z++){
		Ezn1x00[ymax+1][z] = -Ezn1x00[ymax-1][z];
		Ezn1x01[ymax+1][z] = -Ezn1x01[ymax-1][z];
		Ezn1xm0[ymax+1][z] = -Ezn1xm0[ymax-1][z];
		Ezn1xm1[ymax+1][z] = -Ezn1xm1[ymax-1][z];
	}

	// 8回対称
	for(y = 0; y < ymax+1; y++){
		Ezn1x00[y][zmax] = -Ezn1x00[y][zmax-1];
		Ezn1x01[y][zmax] = -Ezn1x01[y][zmax-1];
		Ezn1xm0[y][zmax] = -Ezn1xm0[y][zmax-1];
		Ezn1xm1[y][zmax] = -Ezn1xm1[y][zmax-1];
	}
	for(x = 0; x < xmax+1; x++){
		Ezn1y00[x][zmax] = -Ezn1y00[x][zmax-1];
		Ezn1y01[x][zmax] = -Ezn1y01[x][zmax-1];
		//Ezn1ym0[x][zmax] = -Ezn1ym0[x][zmax-1];
		//Ezn1ym1[x][zmax] = -Ezn1ym1[x][zmax-1];
	}
	for(x = 0; x < xmax; x++){
		Exn1y00[x][zmax+1] = Exn1y00[x][zmax-1];
		Exn1y01[x][zmax+1] = Exn1y01[x][zmax-1];
		//Exn1ym0[x][zmax+1] = Exn1ym0[x][zmax-1];
		//Exn1ym1[x][zmax+1] = Exn1ym1[x][zmax-1];
	}
	for(y = 0; y <= ymax-1; y++){
		Eyn1x00[y][zmax+1] = Eyn1x00[y][zmax-1];
		Eyn1x01[y][zmax+1] = Eyn1x01[y][zmax-1];
		Eyn1xm0[y][zmax+1] = Eyn1xm0[y][zmax-1];
		Eyn1xm1[y][zmax+1] = Eyn1xm1[y][zmax-1];
	}
	/****************************** 対称境界条件 ******************************/




	/****************************** Murの2次の吸収境界条件(Ex) ******************************/

	double u1ax1, u2ax1,u3ax1, u4ax1;
	double u1bx1, u2bx1,u3bx1, u4bx1;
	double u2xa1;

	double velo_dt;

	if(irank != IRANK_MAX){
		for(x = 1; x < xmax; x++){
			for(z = 1; z < zmax+1; z++){
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
			for(z = 1; z < zmax+1; z++){
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

	for(z = 1; z < zmax; z++){

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
			for(z = 1; z < zmax+1; z++){
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
			for(z = 1; z < zmax+1; z++){
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

	for(z = 1; z < zmax+1; z++){
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
		for(z = 1; z < zmax; z++){
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
		for(z = 1; z < zmax; z++){
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
	for(z = 0; z < zmax+1; z++){

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
		for(z = 0; z < zmax+1; z++){
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
			//exn2zm1[x][y] = exn1zm1[x][y];
			//exn1zm1[x][y] = Ex[x][y][zmax-1];
			//exn2zm0[x][y] = exn1zm0[x][y];
			//exn1zm0[x][y] = Ex[x][y][zmax];
		}
	}

	// Ey
	for(x = 0; x < xmax+1; x++){
		for(y = 0; y < ymax; y++){
			Eyn2z00[x][y] = Eyn1z00[x][y];
			Eyn1z00[x][y] = Ey[x][y][0];
			Eyn2z01[x][y] = Eyn1z01[x][y];
			Eyn1z01[x][y] = Ey[x][y][1];
			//eyn2zm1[x][y] = eyn1zm1[x][y];
			//eyn1zm1[x][y] = Ey[x][y][zmax-1];
			//eyn2zm0[x][y] = eyn1zm0[x][y];
			//eyn1zm0[x][y] = Ey[x][y][zmax];
		}
	}
	for(y = 0; y < ymax; y++){
		for(z = 0; z < zmax+1; z++){
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
		for(z = 0; z < zmax; z++){
			Ezn2y00[x][z] = Ezn1y00[x][z];
			Ezn1y00[x][z] = Ez[x][0][z];
			Ezn2y01[x][z] = Ezn1y01[x][z];
			Ezn1y01[x][z] = Ez[x][1][z];
			//ezn2ym1[x][z] = ezn1ym1[x][z];
			//ezn1ym1[x][z] = Ez[x][ymax-1][z];
			//ezn2ym0[x][z] = ezn1ym0[x][z];
			//ezn1ym0[x][z] = Ez[x][ymax][z];
		}
	}
	for(y = 0; y < ymax+1; y++){
		for(z = 0; z < zmax; z++){
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

	//char cell_xy[XMAX][YMAX];
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

	// 全体モデル生成
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


	// YZ平面
	if(irank == intObseInPortNum){ // 入射
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
	}

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


	// ☆ZX平面のモニタ 1/15→1/26これ入れるとプログラムがとまる．
	/*
	if (irank == intObseInPortNum) { // 入射
		y = intObseLenPart1;
		for (z = 0; z < zmax; z++) {
			for (x = 0; x < xmax; x++) {
				cell_zx[z][x] = cell[x][y][z];
			}
		}
		for (z = 0; z < zmax; z++) {
			for (x = 0; x < xmax; x++) {
				fprintf(allmodel_zx, "%d\t", cell_zx[z][x]);
			}
			fprintf(allmodel_zx, "\n");
		}
		fclose(allmodel_zx);
	}*/

	//☆☆ZX平面 //16/1/6
	/*
	if (irank == intObseInPortNum) { // 入射
	x = intObseLenPart1;
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	cell_zx[z][x] = cell[x][y][z];
	}
	}
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	fprintf(allmodel_zx1, "%d\t", cell_zx[y][z]);
	}
	fprintf(allmodel_zx1, "\n");
	}
	fclose(allmodel_zx1);
	}

	if (irank == intObseOutPortNum) { // 出射
	x = intObseLenPart4;
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	cell_zx[z][x] = cell[x][y][z];
	}
	}
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	fprintf(allmodel_zx4, "%d\t", cell_zx[z][x]);
	}
	fprintf(allmodel_zx4, "\n");
	}
	fclose(allmodel_zx4);
	}

	if (irank == intObseCenPortNum) {			// 中央
	x = intObseLenPart7;
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < zmax; x++) {
	cell_zx[z][x] = cell[x][y][z];
	}
	}
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	fprintf(allmodel_zx7, "%d\t", cell_zx[y][z]);
	}
	fprintf(allmodel_zx7, "\n");
	}
	fclose(allmodel_zx7);
	}
	*/




	/****************************** 計算実行時 ******************************/

#else
	/****************************** モデル確認時 ******************************/


	//ZX平面（田村16/01/15）
for (x = 0; x < xmax; x++) {
	for (z = 0; z < zmax; z++) {
		cell_zx[z][x] = cell[x][ymax][z];
	}
}

for (x = 0; x < xmax; x++) {
	for (z = 0; z < zmax; z++) {
		fprintf(model_xz, "%d\t", cell_zx[z][x]);
	}
	fprintf(model_xz, "\n");
}
fclose(model_xz);



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

	//☆☆ZX平面 //16/1/6


	if (irank == intObseInPortNum) { // 入射
	x = intObseLenPart1;
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	cell_zx[z][x] = cell[x][y][z];
	}
	}
	for (z = 0; z < zmax; z++) {
	for (x = 0; x < xmax; x++) {
	fprintf(allmodel_zx1, "%d\t", cell_zx[z][x]);
	}
	fprintf(allmodel_zx1, "\n");
	}
	fclose(allmodel_zx1);
	}

	if(irank == intObseOutPortNum){ // 出射
	x = intObseLenPart4;
	for(z = 0; z < zmax; z++){
	for(x = 0; x < xmax; x++){
	cell_zx[z][x] = cell[x][y][z];
	}
	}
	for(z = 0; z < zmax; z++){
	for(x = 0; x < xmax; x++){
	fprintf(allmodel_zx4, "%d\t", cell_zx[y][z]);
	}
	fprintf(allmodel_zx4, "\n");
	}
	fclose(allmodel_zx4);
	}

	if (irank == intObseCenPortNum) {			// 中央
		x = intObseLenPart7;
		for (z = 0; z < zmax; z++) {
			for (x = 0; x < xmax; x++) {
				cell_zx[z][x] = cell[x][y][z];
			}
		}
		for (z = 0; z < zmax; z++) {
			for (x = 0; x < xmax; x++) {
				fprintf(allmodel_zx7, "%d\t", cell_zx[y][z]);
			}
			fprintf(allmodel_zx7, "\n");
			fclose(allmodel_zx7);
		}

	}
	/****************************** モデル確認時 ******************************/ //☆
#endif

}


void output_field_write(char *dir_name_def) {

	char fname[40], dir_name[50]; 	//ファイル名格納変数
	int node;
	int tag3 = 3;
	int pi1, pj1, pk1;
	MPI_Status status;
	FILE *HZ1;
	//FILE *HZ1, *HZ2;

	pi1 = x_cen;
	pj1 = y_cen;
	pk1 = z_cen;

	if(irank == 0){
		printf("n = %d\n", n);
	}

	for (x = 0; x < xmax; x++) {
		for (y = 0; y < ymax; y++) {
			field_xy[x][y] = Hz[x][y][ex_z_ed - 1]; 	//全てのノードで電磁界成分を2次元配列に格納する．
		}
	}

	// モデル出力ファイルポインタの初期化
	if (irank == IRANK_MIN) {
		sprintf(fname, "/Field_Hz_XY_%d_01.txt", n);
		HZ1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
		for (x = 0; x < xmax; x++) {
			for (y = 0; y < ymax; y++) {
				fprintf(HZ1, "%e\t", field_xy[x][y]);
			}
			fprintf(HZ1, "\n");
		}
	}

	// モデルをホストに送信
	else {
		if (irank != IRANK_MAX) {
			MPI_Send(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 		// ノード0以外のノードがノード0に電磁界成分を送る．
		}
		if (irank == IRANK_MAX) {
			MPI_Send(&field_xy[0][0], (xmax - 1)*(ymax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// ノード0以外のノードがノード0に電磁界成分を送る．
		}
	}

	// 受信したモデルから全モデルを作成
	if (irank == IRANK_MIN) {
		for (node = 1; node < ISIZE; node++) {		// ノード0がノード1から順にデータを受け取り出力していく．
			if (node == IRANK_MAX) {					// ノードisize-1のみ1セル小さく設定しているため条件文で分岐
				MPI_Recv(&field_xy[0][0], (xmax - 1)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for (x = 1; x < xmax - 1; x++) {
					for (y = 0; y < ymax; y++) {
						fprintf(HZ1, "%e\t", field_xy[x][y]);
					}
					fprintf(HZ1, "\n");
				}
			}
			else {
				MPI_Recv(&field_xy[0][0], (xmax)*(ymax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for (x = 1; x < xmax; x++) {
					for (y = 0; y < ymax; y++) {
						fprintf(HZ1, "%e\t", field_xy[x][y]);
					}
					fprintf(HZ1, "\n");
				}
			}
		}

		// ファイルポインタを閉じる
		fclose(HZ1);
	}



	//ZX平面Hy分布（田村16/1/15）☆☆
	FILE *HY1;

	for (x = 0; x < xmax; x++) {
		for (z = 0; z < zmax; z++) {
			field_zx_Hz[z][x] = Hz[x][y_zx][z]; 	//全てのノードで電磁界成分を2次元配列に格納する．
													// y = 0中心，y = 26 シリカと境界(円孔1列目)，y = 89 円孔5列目付近，y = 137 円孔8列目
													//			field_zx[z][x] = Hy[x][y_zx][z];
		}
	}

	/*
	for (x = 0; x < xmax; x++) {//☆☆16/1/22 馬場先生言われたので加えた．
		for (z = 0; z < zmax; z++) {
			field_zx_Hy[z][x] = Hy[x][y_zx][z]; 	//全てのノードで電磁界成分を2次元配列に格納する．
		}
	}*/



	// モデル出力ファイルポインタの初期化
	/*if (irank == IRANK_MIN) {
		sprintf(fname, "/Field_Hz_ZX_%d_01.txt", n);
		HY1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
		for (x = 0; x < xmax; x++) {
			for (z = 0; z < zmax; z++) {
				fprintf(HY1, "%e\t", field_zx_Hz[z][x]);
			}
			fprintf(HY1, "\n");
		}
	}

	// モデルをホストに送信
	else {
		if (irank != IRANK_MAX) {
			MPI_Send(&field_zx_Hz[0][0], (xmax)*(zmax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 		// ノード0以外のノードがノード0に電磁界成分を送る．
		}
		if (irank == IRANK_MAX) {
			MPI_Send(&field_zx_Hz[0][0], (xmax - 1)*(zmax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// ノード0以外のノードがノード0に電磁界成分を送る．
		}
	}*/

	// 受信したモデルから全モデルを作成
/*	if (irank == IRANK_MIN) {
		for (node = 1; node < ISIZE; node++) {		// ノード0がノード1から順にデータを受け取り出力していく．
			if (node == IRANK_MAX) {					// ノードisize-1のみ1セル小さく設定しているため条件文で分岐
				MPI_Recv(&field_zx_Hz[0][0], (xmax - 1)*(zmax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for (x = 1; x < xmax - 1; x++) {
					for (z = 0; z < zmax; z++) {
						fprintf(HY1, "%e\t", field_zx_Hz[z][x]);
					}
					fprintf(HY1, "\n");
				}
			}
			else {
				MPI_Recv(&field_zx_Hz[0][0], (xmax)*(zmax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for (x = 1; x < xmax; x++) {
					for (z = 0; z < zmax; z++) {
						fprintf(HY1, "%e\t", field_zx_Hz[z][x]);
					}
					fprintf(HY1, "\n");
				}
			}
		}

		// ファイルポインタを閉じる
		fclose(HY1);
	}*/


	/*
	// モデル出力ファイルポインタの初期化// Hyの導入に伴い追加
	if (irank == IRANK_MIN) {
		sprintf(fname, "/Field_Hy_ZX_%d_01.txt", n);
		HY1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
		for (x = 0; x < xmax; x++) {
			for (z = 0; z < zmax; z++) {
				fprintf(HY1, "%e\t", field_zx_Hy[z][x]);
			}
			fprintf(HY1, "\n");
		}
	}
	//Field_Hy_ZX
	// モデルをホストに送信
	else {
		if (irank != IRANK_MAX) {
			MPI_Send(&field_zx_Hy[0][0], (xmax)*(zmax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 		// ノード0以外のノードがノード0に電磁界成分を送る．
		}
		if (irank == IRANK_MAX) {
			MPI_Send(&field_zx_Hy[0][0], (xmax - 1)*(zmax), MPI_DOUBLE, 0, tag3, MPI_COMM_WORLD); 	// ノード0以外のノードがノード0に電磁界成分を送る．
		}
	}

	// 受信したモデルから全モデルを作成
	if (irank == IRANK_MIN) {
		for (node = 1; node < ISIZE; node++) {		// ノード0がノード1から順にデータを受け取り出力していく．
			if (node == IRANK_MAX) {					// ノードisize-1のみ1セル小さく設定しているため条件文で分岐
				MPI_Recv(&field_zx_Hy[0][0], (xmax - 1)*(zmax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for (x = 1; x < xmax - 1; x++) {
					for (z = 0; z < zmax; z++) {
						fprintf(HY1, "%e\t", field_zx_Hy[z][x]);
					}
					fprintf(HY1, "\n");
				}
			}
			else {
				MPI_Recv(&field_zx_Hy[0][0], (xmax)*(zmax), MPI_DOUBLE, node, tag3, MPI_COMM_WORLD, &status);
				for (x = 1; x < xmax; x++) {
					for (z = 0; z < zmax; z++) {
						fprintf(HY1, "%e\t", field_zx_Hy[z][x]);
					}
					fprintf(HY1, "\n");
				}
			}
		}

		// ファイルポインタを閉じる
		fclose(HY1);
	}
	*/

	// ☆YZ平面の電界分布の出力
	int x;
	double E_yz;
	FILE *EYZ1, *EYZ2, *EYZ3;
	char fname2[40], fname3[40], fname4[40];

/*	if (irank == intObseInPortNum) { //入射
		x = intObseLenPart1;
		sprintf(fname2, "/Field_E_YZ_%d_01.txt", n);
		EYZ1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname2), "w");

		for (int y = 0; y < ymax; y++) { //矩形導波路断面Y領域判断
			for (int z = 0; z < zmax; z++) {		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				//E_yz = Hz[x][y][y_zx];//→計算するとこれはおかしい

				fprintf(EYZ1, "%e\t", E_yz);
			}
			fprintf(EYZ1, "\n");
		}
		fclose(EYZ1);
	}*/

/*	if (irank == intObseOutPortNum) {			// 出射
		x = intObseLenPart4;
		sprintf(fname3, "/Field_E_YZ_%d_04.txt", n);
		EYZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname3), "w");

		for (int y = 0; y < ymax; y++) { //矩形導波路断面Y領域判断
			for (int z = 0; z < zmax; z++) {		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ2, "%e\t", E_yz);
			}
			fprintf(EYZ2, "\n");
		}
		fclose(EYZ2);
	}

	if (irank == intObseCenPortNum) {			// 出射
		x = intObseLenPart7;
		sprintf(fname4, "/Field_E_YZ_%d_07.txt", n);
		EYZ3 = fopen(strcat(strcpy(dir_name, dir_name_def), fname4), "w");

		for (int y = 0; y < ymax; y++) { //矩形導波路断面Y領域判断
			for (int z = 0; z < zmax; z++) {		// 矩形導波路断面Z領域判断
				E_yz = SQ((Ex[x][y][z] + Ey[x][y][z]));
				fprintf(EYZ3, "%e\t", E_yz);
			}
			fprintf(EYZ3, "\n");
		}
		fclose(EYZ3);
	}*/

	// 2013/01/22 最終ステップでのYZ平面のHz成分を出力
	//if (n == Nmax){
	//	if (irank == intObseInPortNum){
	//		sprintf(fname, "/Field_Hz_YZ_%d_03.txt", n);
	//		HZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
	//		for(y = 0; y < ymax; y++){
	//			for(z = 0; z < zmax; z++){
	//				fprintf(HZ2, "%e\t", Hz[intObseLenPart3][y][z]);
	//			}
	//			fprintf(HZ2, "\n");
	//		}
	//		fclose(HZ2);
	//	}
	//	if (irank == intObseOutPortNum) {
	//		sprintf(fname, "/Field_Hz_YZ_%d_06.txt", n);
	//		HZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
	//		for(y = 0; y < ymax; y++){
	//			for(z = 0; z < zmax; z++){
	//				fprintf(HZ2, "%e\t", Hz[intObseLenPart6][y][z]);
	//			}
	//			fprintf(HZ2, "\n");
	//		}
	//		fclose(HZ2);
	//	}
	//}



//☆☆☆ZX平面の電界分布の出力16/1/5 プログラム作成開始 ここで別のモニタとして機能

/*	double E_zx;
	char fname10[40];
	FILE *EZX1, *EZX2, *EZX3;


	char fname5[40], fname6[40], fname7[40];
	if (irank == intObseInPortNum) { //入射
	x = intObseLenPart1;
	sprintf(fname5, "/Field_E_ZX_%d_01.txt", n);
	EZX1 = fopen(strcat(strcpy(dir_name, dir_name_def), fname5), "w");

	for (int z = 0; z < zmax; z++) { //矩形導波路断面Z領域判断
	for (int x = 0; x < xmax; x++) {		// 矩形導波路断面X領域判断
	E_zx = SQ((Ey[x][y][z] + Ez[x][y][z]));
	fprintf(EZX1, "%e\t", E_zx);
	}
	fprintf(EZX1, "\n");
	}
	fclose(EZX1);
	}


		if(irank == intObseOutPortNum){			// 出射
			x = intObseLenPart4;
			sprintf(fname6, "/Field_E_ZX_%d_04.txt", n);
			EZX2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname6), "w");

	for (int z = 0; z < zmax; z++) { //矩形導波路断面Z領域判断
	for (int x = 0; x < xmax; x++) {		// 矩形導波路断面X領域判断
	E_zx = SQ((Ey[x][y][z] + Ez[x][y][z]));
		fprintf(EZX2, "%e\t", E_zx);
		}
		fprintf(EZX2, "\n");
		}
		fclose(EZX2);
	}



	if (irank == intObseCenPortNum) {			// 出射//ここだけ使えばいい
	x = intObseLenPart7;
	sprintf(fname7, "/Field_E_ZX_%d_07.txt", n);
	EZX3 = fopen(strcat(strcpy(dir_name, dir_name_def), fname7), "w");

	for (int z = 0; z < zmax; z++) { //矩形導波路断面Z領域判断
	for (int x = 0; x < xmax; x++) {		// 矩形導波路断面X領域判断
	E_zx = SQ((Ey[x][y][z] + Ez[x][y][z]));
	fprintf(EZX3, "%e\t", E_zx);
	}
	fprintf(EZX3, "\n");
	}
	fclose(EZX3);
	}

*/

	//// 2013/01/22 最終ステップでのYZ平面のHz成分を出力
	//if (n == Nmax){
	//	if (irank == intObseInPortNum){
	//		sprintf(fname, "/Field_Hz_YZ_%d_03.txt", n);
	//		HZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
	//		for(y = 0; y < ymax; y++){
	//			for(z = 0; z < zmax; z++){
	//				fprintf(HZ2, "%e\t", Hz[intObseLenPart3][y][z]);
	//			}
	//			fprintf(HZ2, "\n");
	//		}
	//		fclose(HZ2);
	//	}
	//	if (irank == intObseOutPortNum) {
	//		sprintf(fname, "/Field_Hz_YZ_%d_06.txt", n);
	//		HZ2 = fopen(strcat(strcpy(dir_name, dir_name_def), fname), "w");
	//		for(y = 0; y < ymax; y++){
	//			for(z = 0; z < zmax; z++){
	//				fprintf(HZ2, "%e\t", Hz[intObseLenPart6][y][z]);
	//			}
	//			fprintf(HZ2, "\n");
	//		}
	//		fclose(HZ2);
	//	}
	//}
}

//☆ファイル出力
void output_field(char *dir_name_def){

	//double field_xy[XMAX][YMAX]; 	// Hz-field のファイル出力 (面垂直方向の磁界成分)

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


/*//☆☆☆ファイル出力
void output_field(char *dir_name_def) {

// Hz-field のファイル出力 (面導波路垂直方向の磁界成分)16/1/5 プログラム作成開始
// もしかしたらこれは違うか？

if (n <= Nmax - Fcut) {
// 動作確認のためのファイル出力
if (n == Ncheck) {
output_field_write(dir_name_def);
}

// 定期的なファイル出力
if (n % Ncutfield == 0) {
output_field_write(dir_name_def);
}
}
if ((n >= Nmax - Fcut) && (n <= Nmax)) {

// 安定点でのファイル出力
if (n % Ncutfield2 == 0) {
output_field_write(dir_name_def);
}
}
}*/


void calc_poynting_powerHz(){				// ☆透過スペクトル計算用の観測面中央でのHzの出力

	//入力部分での poynting power (x方向)
	//double Eyin, Ezin, Hyin, Hzin; 		// 各成分の積算保存変数
	//double EyHz, EzHy; 					// 積分値を掛け算を保存する変数とセル中央の値

	double pmax_01 = 0;		// 入力パワーの最大値を記録する変数
	double pmax_03 = 0;		// 出力パワーの最大値を記録
	double pmin_01 = 0;		// 出力パワーの最小値を記録
	double pmin_03 = 0;		// 出力パワーの最小値を記録

	int x = 0;

	if(irank == intObseInPortNum){ //入射
		x = intObseLenPart1;
		y = ymax - 1;
		z = zmax - 1;

		fprintf(fpHz1, "%e\n", Hz[x][y][z]);

	}

	if(irank == intObseOutPortNum){			// 出射
		x = intObseLenPart4;
		y = ymax - 1;
		z = zmax - 1;

		fprintf(fpHz5, "%e\n", Hz[x][y][z]);

	}
}


/*
void calc_poynting_power(){				// ☆吉田作成パワー評価プログラム

	//入力部分での poynting power (x方向)
	double Eyin, Ezin, Hyin, Hzin; 		// 各成分の積算保存変数
	double EyHz, EzHy; 					// 積分値の掛け算を保存する変数とセル中央の値

	double pmax_01 = 0;		// 入力パワーの最大値を記録する変数
	double pmax_03 = 0;		// 出力パワーの最大値を記録
	double pmin_01 = 0;		// 出力パワーの最小値を記録 //入射の間違えでは？
	double pmin_03 = 0;		// 出力パワーの最小値を記録

	int x = 0;

	if(irank == intObseInPortNum){ //入射

		if(n == Nmax - Tcut){ //ある時刻までの入出力パワーの最大値と最小値を保持する変数．16/2/1 ある時刻はいつまでか？
			powermax_in = 0.0;
			powermin_in = 0.0;
			powermax_out = 0.0;
			powermin_out = 0.0;
		}

		for(x = intObseLenPart1; x < intObseLenPart2; x++){

			// 初期化
			Eyin = 0; 	Ezin = 0; 	Hyin = 0; 	Hzin = 0; 	EyHz = 0; 	EzHy = 0;

			// 観測面の修正(2013/8/8) /
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
					//	for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
					// 観測面の修正(2013/8/8)


					Eyin = 0.25 * (Ey[x][y][z] + Ey[x+1][y][z] + Ey[x][y][z+1] + Ey[x+1][y][z+1]);
					Hzin = 0.50 * (Hz[x][y][z] + Hz[x][y][z+1]);
					Ezin = 0.25 * (Ez[x][y][z] + Ez[x+1][y][z] + Ez[x][y+1][z] + Ez[x+1][y+1][z]);
					Hyin = 0.50 * (Hy[x][y][z] + Hy[x][y+1][z]);

					EyHz += Eyin * Hzin;
					EzHy += Ezin * Hyin;
				}
			}

			//ポインティングパワーをファイルに保存
			if(x == intObseLenPart3){
				fprintf(fppoynt1, "%e\n", EyHz - EzHy);
			}
			//if(x == intObseLenPartHz1){
			//	fprintf(fppoynt1h, "%e\n", EyHz - EzHy);
			//}

			if(n >= Nmax - Tcut){
				if(pmax_01 < EyHz - EzHy){
					pmax_01 = EyHz - EzHy;
				}

				if(pmin_01 > EyHz - EzHy){
					pmin_01 = EyHz - EzHy;
				}
			}
		}


		if(n >= Nmax - Tcut){
			if(powermax_in < pmax_01){
				powermax_in = pmax_01;
			}
			if(powermin_in > pmin_01){
				powermin_in = pmin_01;
			}
		}
		if(n == Nmax){
			fprintf(avpoynt1, "%s", "Input_Power\t");
			fprintf(avpoynt1, "%e\n", powermax_in);
			fprintf(avpoynt1, "%s", "Reflection\t");
			fprintf(avpoynt1, "%e\n", powermin_in);
			fprintf(avpoynt1, "%s", "Degree_of_Reflection\t");
			fprintf(avpoynt1, "%e\n", powermin_in/powermax_in);
		}
	}

	if(irank == intObseOutPortNum){			// 出射

		if(n == Nmax - Tcut){ //ある時刻までの入出力パワーの最大値と最小値を保持する変数
			powermax_in = 0.0;
			powermin_in = 0.0;
			powermax_out = 0.0;
			powermin_out = 0.0;
		}

		for(x = intObseLenPart4; x < intObseLenPart5; x++){//ここでエラー発生()16/1/6午前

			//初期化
			Eyin = 0; 	Ezin = 0; 	Hyin = 0; 	Hzin = 0; 	EyHz = 0; 	EzHy = 0;

			// 観測面の修正(2013/8/8) /
			for(int y = ymax - intObseWid; y < ymax; y++){ // 矩形導波路断面Y領域判断．
				for(int z = zmax - intObseHeig; z < zmax; z++){		//矩形導波路断面Z領域判断．
					//for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
					//	for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
					//観測面の修正(2013/8/8) /

					Eyin = 0.25 * (Ey[x][y][z] + Ey[x+1][y][z] + Ey[x][y][z+1] + Ey[x+1][y][z+1]);
					Hzin = 0.50 * (Hz[x][y][z] + Hz[x][y][z+1]);
					Ezin = 0.25 * (Ez[x][y][z] + Ez[x+1][y][z] + Ez[x][y+1][z] + Ez[x+1][y+1][z]);
					Hyin = 0.50 * (Hy[x][y][z] + Hy[x][y+1][z]);

					//if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
					//	cell[x][y][z] = 4; 					// 中央点確認用
					//}

					EyHz += Eyin * Hzin;
					EzHy += Ezin * Hyin;					// ポインティングパワーを各セル毎に足す
				}
			}

			if(x == intObseLenPart6){
				fprintf(fppoynt5, "%e\n", EyHz - EzHy); 		// ポインティングパワーをファイルに保存
			}
			//if(x == intObseLenPartHz5){
			//	fprintf(fppoynt5h, "%e\n", EyHz - EzHy);
			//}

			if(n >= Nmax - Tcut){
				if(pmax_03 < EyHz - EzHy){
					pmax_03 = EyHz - EzHy;
				}
				if(pmin_03 > EyHz - EzHy){
					pmin_03 = EyHz - EzHy;
				}
			}
		}
		if(n >= Nmax - Tcut){
			if(powermax_out < pmax_03){
				powermax_out = pmax_03;
			}

			if(powermin_out > pmin_03){
				powermin_out = pmin_03;
			}
		}
		if(n == Nmax){
			fprintf(avpoynt5, "%s", "Transmission\t");
			fprintf(avpoynt5, "%e\n", powermax_out);
			fprintf(avpoynt5, "%s", "Reflection\t");
			fprintf(avpoynt5, "%e\n\n", powermin_out);
		}
	}
}*/


void calc_poynting_power() {				// ☆吉田作成パワー評価プログラム改 途中経過(not～Nmax)測定用 16/2/1

											//入力部分での poynting power (x方向)
	double Eyin, Ezin, Hyin, Hzin; 		// 各成分の積算保存変数
	double EyHz, EzHy; 					// 積分値の掛け算を保存する変数とセル中央の値

	double pmax_01 = 0;		// 入力パワーの最大値を記録する変数
	double pmax_03 = 0;		// 出力パワーの最大値を記録
	double pmin_01 = 0;		// 出力パワーの最小値を記録 //入射の間違えでは？
	double pmin_03 = 0;		// 出力パワーの最小値を記録

	int x = 0;

	if (irank == intObseInPortNum) { //入射

		if (n == Nmax - Tcut) { //ある時刻までの入出力パワーの最大値と最小値を保持する変数．16/2/1 ある時刻はいつまでか？
			powermax_in = 0.0;
			powermin_in = 0.0;
			powermax_out = 0.0;
			powermin_out = 0.0;
		}

		for (x = intObseLenPart1; x < intObseLenPart2; x++) {

			// 初期化
			Eyin = 0; 	Ezin = 0; 	Hyin = 0; 	Hzin = 0; 	EyHz = 0; 	EzHy = 0;

			/****************************** 観測面の修正(2013/8/8) ******************************/
			for (int y = ymax - intObseWid; y < ymax; y++) { // 矩形導波路断面Y領域判断．
				for (int z = zmax - intObseHeig; z < zmax; z++) {		//矩形導波路断面Z領域判断．
																		//for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
																		//	for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
																		/****************************** 観測面の修正(2013/8/8) ******************************/


					Eyin = 0.25 * (Ey[x][y][z] + Ey[x + 1][y][z] + Ey[x][y][z + 1] + Ey[x + 1][y][z + 1]);
					Hzin = 0.50 * (Hz[x][y][z] + Hz[x][y][z + 1]);
					Ezin = 0.25 * (Ez[x][y][z] + Ez[x + 1][y][z] + Ez[x][y + 1][z] + Ez[x + 1][y + 1][z]);
					Hyin = 0.50 * (Hy[x][y][z] + Hy[x][y + 1][z]);

					EyHz += Eyin * Hzin;
					EzHy += Ezin * Hyin;
				}
			}

			//ポインティングパワーをファイルに保存
			if (x == intObseLenPart3) {
				fprintf(fppoynt1, "%e\n", EyHz - EzHy);
			}
			//if(x == intObseLenPartHz1){
			//	fprintf(fppoynt1h, "%e\n", EyHz - EzHy);
			//}

			if (n >= Nmax - Tcut) {
				if (pmax_01 < EyHz - EzHy) {
					pmax_01 = EyHz - EzHy;
				}

				if (pmin_01 > EyHz - EzHy) {
					pmin_01 = EyHz - EzHy;
				}
			}
		}


		if (n >= Nmax - Tcut) {
			if (powermax_in < pmax_01) {
				powermax_in = pmax_01;
			}
			if (powermin_in > pmin_01) {
				powermin_in = pmin_01;
			}
		}
		if (n == Nmax) {
			fprintf(avpoynt1, "%s", "Input_Power\t");
			fprintf(avpoynt1, "%e\n", powermax_in);
			fprintf(avpoynt1, "%s", "Reflection\t");
			fprintf(avpoynt1, "%e\n", powermin_in);
			fprintf(avpoynt1, "%s", "Degree_of_Reflection\t");
			fprintf(avpoynt1, "%e\n", powermin_in / powermax_in);
		}
	}

	if (irank == intObseOutPortNum) {			// 出射

		if (n == Nmax - Tcut) { //ある時刻までの入出力パワーの最大値と最小値を保持する変数
			powermax_in = 0.0;
			powermin_in = 0.0;
			powermax_out = 0.0;
			powermin_out = 0.0;
		}

		for (x = intObseLenPart4; x < intObseLenPart5; x++) {//ここでエラー発生()16/1/6午前

															 //初期化
			Eyin = 0; 	Ezin = 0; 	Hyin = 0; 	Hzin = 0; 	EyHz = 0; 	EzHy = 0;

			/****************************** 観測面の修正(2013/8/8) ******************************/
			for (int y = ymax - intObseWid; y < ymax; y++) { // 矩形導波路断面Y領域判断．
				for (int z = zmax - intObseHeig; z < zmax; z++) {		//矩形導波路断面Z領域判断．
																		//for(int y = 0; y < ymax; y++){ //矩形導波路断面Y領域判断
																		//	for(int z = 0; z < zmax; z++){		// 矩形導波路断面Z領域判断
																		/****************************** 観測面の修正(2013/8/8) ******************************/

					Eyin = 0.25 * (Ey[x][y][z] + Ey[x + 1][y][z] + Ey[x][y][z + 1] + Ey[x + 1][y][z + 1]);
					Hzin = 0.50 * (Hz[x][y][z] + Hz[x][y][z + 1]);
					Ezin = 0.25 * (Ez[x][y][z] + Ez[x + 1][y][z] + Ez[x][y + 1][z] + Ez[x + 1][y + 1][z]);
					Hyin = 0.50 * (Hy[x][y][z] + Hy[x][y + 1][z]);

					//if((y == YMAX-1) && (z == (intSlabCen-1))){		// 活性層断面中央点の時を記録
					//	cell[x][y][z] = 4; 					// 中央点確認用
					//}

					EyHz += Eyin * Hzin;
					EzHy += Ezin * Hyin;					// ポインティングパワーを各セル毎に足す
				}
			}

			if (x == intObseLenPart6) {
				fprintf(fppoynt5, "%e\n", EyHz - EzHy); 		// ポインティングパワーをファイルに保存
			}
			//if(x == intObseLenPartHz5){
			//	fprintf(fppoynt5h, "%e\n", EyHz - EzHy);
			//}

			if (n >= Nmax - Tcut) {
				if (pmax_03 < EyHz - EzHy) {
					pmax_03 = EyHz - EzHy;
				}
				if (pmin_03 > EyHz - EzHy) {
					pmin_03 = EyHz - EzHy;
				}
			}
		}
		if (n >= Nmax - Tcut) {
			if (powermax_out < pmax_03) {
				powermax_out = pmax_03;
			}

			if (powermin_out > pmin_03) {
				powermin_out = pmin_03;
			}
		}
		if (n == Nmax) {
			fprintf(avpoynt5, "%s", "Transmission\t");
			fprintf(avpoynt5, "%e\n", powermax_out);
			fprintf(avpoynt5, "%s", "Reflection\t");
			fprintf(avpoynt5, "%e\n\n", powermin_out);
		}
	}
}





void mcircle(int x_circ, int y_circ, int z_circ, int type) {//☆☆ ここで円孔の大きさ制御か？

	double R;

	//半径セル数の計算
	if (type == 1)	R = ((dblRadius*1.0e10) / (dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています
	else if (type == 2)	R = ((dblRadius2*1.0e10) / (dx*1.0e10));
	else if (type == 3)	R = ((dblRadius3*1.0e10) / (dx*1.0e10));
	else if (type == 5)	R = ((dblRadius5*1.0e10) / (dx*1.0e10));
	else if (type == 6)	R = ((dblRadius6*1.0e10) / (dx*1.0e10));
	else if (type == 7)	R = ((dblRadius7*1.0e10) / (dx*1.0e10));
	else if (type == 8)	R = ((dblRadius8*1.0e10) / (dx*1.0e10));

	else if (type == 9)	R = ((dblRadius9*1.0e10) / (dx*1.0e10));
	else if (type == 10)	R = ((dblRadius10*1.0e10) / (dx*1.0e10));
	else if (type == 11)	R = ((dblRadius11*1.0e10) / (dx*1.0e10));
	else if (type == 12)	R = ((dblRadius12*1.0e10) / (dx*1.0e10));
	else if (type == 13)	R = ((dblRadius13*1.0e10) / (dx*1.0e10));
	else if (type == 14)	R = ((dblRadius14*1.0e10) / (dx*1.0e10));
	else if (type == 15)	R = ((dblRadius15*1.0e10) / (dx*1.0e10));
	else			R = ((dblRadius4*1.0e10) / (dx*1.0e10));

	rightquartercircle1(x_circ, y_circ, z_circ, type, R);
	leftquartercircle1(x_circ - 1, y_circ, z_circ, type, R);
	rightquartercircle2(x_circ, y_circ - 1, z_circ, type, R);
	leftquartercircle2(x_circ - 1, y_circ - 1, z_circ, type, R);
}


void halfcircle(int x_circ, int y_circ, int z_circ, int type) {//☆☆多分これは使っていない

	double R;

	//半径セル数の計算
	if (type == 1)	R = ((dblRadius*1.0e10) / (dx*1.0e10)); 		//計算誤差を防ぐために桁上げしています
	else if (type == 2)	R = ((dblRadius2*1.0e10) / (dx*1.0e10));
	else if (type == 3)	R = ((dblRadius3*1.0e10) / (dx*1.0e10));
	else if (type == 5)	R = ((dblRadius5*1.0e10) / (dx*1.0e10));
	else if (type == 6)	R = ((dblRadius6*1.0e10) / (dx*1.0e10));
	else if (type == 7)	R = ((dblRadius7*1.0e10) / (dx*1.0e10));
	else if (type == 8)	R = ((dblRadius8*1.0e10) / (dx*1.0e10));

	else if (type == 9)	R = ((dblRadius9*1.0e10) / (dx*1.0e10));
	else if (type == 10)	R = ((dblRadius10*1.0e10) / (dx*1.0e10));
	else if (type == 11)	R = ((dblRadius11*1.0e10) / (dx*1.0e10));
	else if (type == 12)	R = ((dblRadius12*1.0e10) / (dx*1.0e10));
	else if (type == 13)	R = ((dblRadius13*1.0e10) / (dx*1.0e10));
	else if (type == 14)	R = ((dblRadius14*1.0e10) / (dx*1.0e10));
	else if (type == 15)	R = ((dblRadius15*1.0e10) / (dx*1.0e10));
	else			R = ((dblRadius4*1.0e10) / (dx*1.0e10));

	rightquartercircle2(x_circ, y_circ - 1, z_circ, type, R);
	leftquartercircle2(x_circ - 1, y_circ - 1, z_circ, type, R);
}

void rightquartercircle1(int x_circ, int y_circ, int z_circ, int type, double R) {

	int x, y, Ie, Je;
	double r;

	Ie = (int)(x_circ + R/*-1*/);
	Je = (int)(y_circ + R/*-1*/);
	for (x = x_circ; x <= Ie; x++) {
		for (y = y_circ; y <= Je; y++) {
			r = sqrt(double((x - x_circ + 1) * (x - x_circ + 1) + (y - y_circ + 1) * (y - y_circ + 1))) - 0.5;
			if (r <= R) {
				if (type == 1 || type == 3 || type == 5 || type == 6 || type == 7 || type == 8 || type == 9 || type == 10 || type == 11 || type == 12 || type == 13 || type == 14 || type == 15) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else if (type == 2) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}


void leftquartercircle1(int x_circ, int y_circ, int z_circ, int type, double R) {

	int x, y, Ie, Je;
	double r;

	Ie = (int)(x_circ - R/*+1*/);
	Je = (int)(y_circ + R/*-1*/);
	for (x = x_circ; x >= Ie; x--) {
		for (y = y_circ; y <= Je; y++) {
			r = sqrt(double((x - x_circ - 1) * (x - x_circ - 1) + (y - y_circ + 1) * (y - y_circ + 1))) - 0.5;
			if (r <= R) {
				if (type == 1 || type == 3 || type == 5 || type == 6 || type == 7 || type == 8 || type == 9 || type == 10 || type == 11 || type == 12 || type == 13 || type == 14 || type == 15) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else if (type == 2) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}


void rightquartercircle2(int x_circ, int y_circ, int z_circ, int type, double R) {

	int x, y, Ie, Je;
	double r;

	Ie = (int)(x_circ + R/*-1*/);
	Je = (int)(y_circ - R/*+1*/);
	for (x = x_circ; x <= Ie; x++) {
		for (y = y_circ; y >= Je; y--) {
			r = sqrt(double((x - x_circ + 1) * (x - x_circ + 1) + (y - y_circ - 1) * (y - y_circ - 1))) - 0.5;
			if (r <= R) {
				if (type == 1 || type == 3 || type == 5 || type == 6 || type == 7 || type == 8 || type == 9 || type == 10 || type == 11 || type == 12 || type == 13 || type == 14 || type == 15) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else if (type == 2) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}
}

void leftquartercircle2(int x_circ, int y_circ, int z_circ, int type, double R) {

	int x, y, Ie, Je;
	double r;

	Ie = (int)(x_circ - R/*+1*/);
	Je = (int)(y_circ - R/*+1*/);
	for (x = x_circ; x >= Ie; x--) {
		for (y = y_circ; y >= Je; y--) {
			r = sqrt(double((x - x_circ - 1) * (x - x_circ - 1) + (y - y_circ - 1) * (y - y_circ - 1))) - 0.5;
			if (r <= R) {
				if (type == 1 || type == 3 || type == 5 || type == 6 || type == 7 || type == 8 || type == 9 || type == 10 || type == 11 || type == 12 || type == 13 || type == 14 || type == 15) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else if (type == 2) {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX2;
					ALL_epsilonx[x][y][z_circ] = epsilon2;
					ALL_epsilony[x][y][z_circ] = epsilon2;
					ALL_epsilonz[x][y][z_circ] = epsilon2;
				}
				else {
					ALL_cell[x][y][z_circ] = CIRCLE_REF_INDEX3;
				}
			}
		}
	}



}
