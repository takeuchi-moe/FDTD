/*****************************************************************************/
// 一つの計算機での解析領域
/*****************************************************************************/
#define NUM_EX 3
#define XMAX (INT_DIV(XMAX_ALL, 3)*NUM_EX + 1) //(INT_DIV(XMAX_ALL, NODE) + 1)
#define YMAX YMAX_ALL
#define ZMAX ZMAX_ALL
#define X_LSPCW (INT_DIV(XMAX_ALL, 3) + 1)

/*****************************************************************************/
// 励振関数
/*****************************************************************************/
//char *dir_name[] = {"1550", "1580"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1550"}; 		// 励振関数の波長 [nm]
//char *dir_name[] = {"1555", "1560", "1570", "1575"}; 		// 励振関数の波長 [nm]
double lambda; 											// 励振関数の波長 [m] (プログラム中で dir_name をdouble化して代入)
double omega0; 											// 励振関数の角周波数 [s^-1]

// ガウシアンパルス
double sigma; 											// 広がり幅を決定する定数
//static const double delta_omega = 0.1; 						// 中心周波数で規格した半値全幅
//static const int Npeak = 3000; 								// ピークステップ数


/*****************************************************************************/
// Model掃き出し用の定数
/*****************************************************************************/
#define CLAD 0									// クラッド
#define CORE 1									// コア
#define GaInAsP	2								// GaInAsP
#define AIR_GaInAsP	3							// CLAD/GaInAsP
#define EXITATION 30							// 励振点
#define OBSERVATION 20							// 観測点
#define CIRCLE_REF_INDEX	CLAD				//関数mcircleで書き込む数字
#define CIRCLE_REF_INDEX2	CLAD		//関数mcircle2で書き込む数字
#define CIRCLE_REF_INDEX3	2				//場所確認用ドットの色指定


/*****************************************************************************/
// フォトニック結晶モデル用宣言
/*****************************************************************************/
//厚さ方向パラメータ
struct PNUM {int X; int Y; }; 					//円柱の中心座標を与える変数.sankakuで使用


/*****************************************************************************/
// グローバル変数 (MPIの通信ではグローバル変数を使用しないとエラーが生じる模様)
/*****************************************************************************/

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION

// 電磁界 (XMAX+1は並列計算での"のりしろ"が必要なため，YMAX+1, ZMAX+1は対称境界条件で必要なため)
double Ex[XMAX+1][YMAX+1][ZMAX+1];
double Ey[XMAX+1][YMAX+1][ZMAX+1];
double Ez[XMAX+1][YMAX+1][ZMAX+1];
double Hx[XMAX+1][YMAX+1][ZMAX+1];
double Hy[XMAX+1][YMAX+1][ZMAX+1];
double Hz[XMAX+1][YMAX+1][ZMAX+1];

// 誘電率 (YMAX+1, ZMAX+1は対称境界条件で必要なため)
double ALL_epsilonx[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];
double ALL_epsilony[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];
double ALL_epsilonz[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];
double epsilonx[XMAX][YMAX+1][ZMAX+1];
double epsilony[XMAX][YMAX+1][ZMAX+1];
double epsilonz[XMAX][YMAX+1][ZMAX+1];
int ALL_cell[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];		// 実際はALL_cell[XMAX_ALL+1][YMAX_ALL][ZMAX_ALL]で良いが，プログラムの都合上誘電体の配列と同じ次元にした
int cell[XMAX][YMAX+1][ZMAX+1];							// 実際はcell[XMAX][YMAX][ZMAX]で良いが，プログラムの都合上誘電体の配列と同じ次元にした
double epsilon_xy[XMAX+1][YMAX+1], epsilon_yz[YMAX+1][ZMAX+1], epsilon_zx[XMAX+1][ZMAX+1];
double epsilon_zx2[XMAX+1][ZMAX+1];
int cell_xy[XMAX][YMAX], cell_yz[YMAX][ZMAX];

/*吸収境界条件適用のときに使う電界の配列*/
double Exn2y00[XMAX+1][ZMAX+1], 	Exn1y00[XMAX+1][ZMAX+1+1], 	Exn2y01[XMAX+1][ZMAX+1], 	Exn1y01[XMAX+1][ZMAX+1+1];
//double exn2ym1[XMAX+1][ZMAX+1], 	exn1ym1[XMAX+1][ZMAX+1+1], 	exn2ym0[XMAX+1][ZMAX+1], 	exn1ym0[XMAX+1][ZMAX+1+1];
double Exn2z00[XMAX+1][YMAX+1], 	Exn1z00[XMAX+1][YMAX+1+1], 	Exn2z01[XMAX+1][YMAX+1], 	Exn1z01[XMAX+1][YMAX+1+1];
//double exn2zm1[XMAX+1][YMAX+1], 	exn1zm1[XMAX+1][YMAX+1+1], 	exn2zm0[XMAX+1][YMAX+1], 	exn1zm0[XMAX+1][YMAX+1+1];

double Eyn2z00[XMAX+1][YMAX], 	Eyn1z00[XMAX+1][YMAX+1], 	Eyn2z01[XMAX+1][YMAX], 	Eyn1z01[XMAX+1][YMAX+1];
//double eyn2zm1[XMAX+1][YMAX], 	eyn1zm1[XMAX+1][YMAX+1], 	eyn2zm0[XMAX+1][YMAX], 	eyn1zm0[XMAX+1][YMAX+1];
double Eyn2x00[YMAX][ZMAX+1], 	Eyn1x00[YMAX+1][ZMAX+1+1], 	Eyn2x01[YMAX][ZMAX+1], 	Eyn1x01[YMAX+1][ZMAX+1+1];
double Eyn2xm1[YMAX][ZMAX+1], 	Eyn1xm1[YMAX+1][ZMAX+1+1], 	Eyn2xm0[YMAX][ZMAX+1], 	Eyn1xm0[YMAX+1][ZMAX+1+1];

double Ezn2y00[XMAX+1][ZMAX], 	Ezn1y00[XMAX+1][ZMAX+1], 	Ezn2y01[XMAX+1][ZMAX], 	Ezn1y01[XMAX+1][ZMAX+1];
//double ezn2ym1[XMAX+1][ZMAX], 	ezn1ym1[XMAX+1][ZMAX+1], 	ezn2ym0[XMAX+1][ZMAX], 	ezn1ym0[XMAX+1][ZMAX+1];
double Ezn2x00[YMAX+1][ZMAX], 	Ezn1x00[YMAX+1+1][ZMAX+1], 	Ezn2x01[YMAX+1][ZMAX], 	Ezn1x01[YMAX+1+1][ZMAX+1];
double Ezn2xm1[YMAX+1][ZMAX], 	Ezn1xm1[YMAX+1+1][ZMAX+1], 	Ezn2xm0[YMAX+1][ZMAX], 	Ezn1xm0[YMAX+1+1][ZMAX+1];

// 磁界分布出力用
double field_xy[XMAX][YMAX]; 	// Hz-field のファイル出力 (面垂直方向の磁界成分)
double field_xz[XMAX][ZMAX]; 	// Hz-field のファイル出力 (面垂直方向の磁界成分)
double field_yz[YMAX][ZMAX];	// Hz-field のファイル出力 (面方向の磁界成分)

int x, y, z; 				// 離散座標
int xmax, ymax, zmax; 			// 空間分割の最大値
int xmax_all, ymax_all, zmax_all; //分割前の最大値
//int n, Nmax; 				//n: 時間ステップ，Nmax: 時間ステップの最大値
int n; 					//時間ステップ数
//int icut, jcut, kcut; 		//icut, jcut, kcut: フィールド出力のときに使う
//int PrintStat;
//int PrintEnd;
//int source_i1, source_i2;
//int source_j1, source_j2;
//int source_k;
int x_cen, y_cen, z_cen; 			//x_cen, y_cen, z_cen: 解析空間の中心離散座標
int x_model_cen, y_model_cen; 	//x_model_cen, y_model_cen:モデルの中心離散座標
//char fname[30];
int irank, isize;

// 電磁界計算に使う定数の宣言
static const double cnstHxyz = dt / MU0; 					//磁界の計算に使う定数		cnstHxyz: Magnetic field calculation


///*励振用定数の宣言*/
//double omega0; 	//omega0: 励振関数の角周波数
//double sigma; 	//ガウシアンパルスのための広がり幅を決定する定数
//int Npeak; 		//ガウシアンパルスのピークステップ数

/*モデル設定に使う定数の宣言 -- 関数modeling()以外でも使う定数があるのでここで宣言する．*/
double b; //ディスク厚さ
int bk; //ディスク厚さの離散値

double w; 		//支柱幅の半分
int wij; 		//支柱幅の半分の離散値
double h; 		//支柱の高さ
int hk; 			//支柱の高さの離散値
//double np; 		//支柱の屈折率

//入出力パワーの最大値と最小値を記録する変数
double powermax_in;
double powermin_in;
double powermax_out;
double powermin_out;

double topy1, topy5, tohz1, tohz5;
double topy1h, topy5h, tohz1h, tohz5h;

void mcircle(int, int , int, int);  //make circle function
void halfcircle(int, int , int, int); 							//make circle function
void rightquartercircle1(int, int, int, int, double);
void leftquartercircle1(int, int, int, int, double);
void rightquartercircle2(int, int, int, int, double);
void leftquartercircle2(int, int, int, int, double);

/////////////////////////////////////////////////////
/*観測点*/
int inputi, inputj, inputk;
int outputi, outputj, outputk;
double pyn_in, pyn_out;


/*ファイルポインタ宣言*/
//モデル			//Ez				//Hz
FILE *model_xy; 		FILE *fpfez_xy; 		FILE *fpfhz_xy;
FILE *model_yz; 		FILE *fpfez_yz; 		FILE *fpfhz_yz;
FILE *model_xz; 		FILE *fpfez_zx; 		FILE *fpfhz_zx;
FILE *model_xy2; 	FILE *fpfez2_xy; 	FILE *fpfhz2_xy;
FILE *allmodel_xy;
FILE *allmodel_yz1, *allmodel_yz4, *allmodel_yz7;

FILE *model_xy_1, *model_xy_2;

//誘電率
FILE *fpepsilonx;
FILE *fpallepsilonx;
FILE *fpepsilony;
FILE *fpepsilonz;
FILE *fpepsilony2;
FILE *fpepsilonz2;
FILE *fpAllEpsilon;
FILE *fpEpsilon;


FILE *fpex; 		FILE *fphx;
FILE *fpey; 		FILE *fphy;
FILE *fpez; 		FILE *fphz;

FILE *fpparameter; 	//計算パラメータ保存ファイル


#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
