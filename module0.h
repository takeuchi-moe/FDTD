/*****************************************************************************/
// ��̌v�Z�@�ł̉�͗̈�
/*****************************************************************************/
#define NUM_EX 3
#define XMAX (INT_DIV(XMAX_ALL, 3)*NUM_EX + 1) //(INT_DIV(XMAX_ALL, NODE) + 1)
#define YMAX YMAX_ALL
#define ZMAX ZMAX_ALL
#define X_LSPCW (INT_DIV(XMAX_ALL, 3) + 1)

/*****************************************************************************/
// ��U�֐�
/*****************************************************************************/
//char *dir_name[] = {"1550", "1580"}; 		// ��U�֐��̔g�� [nm]
//char *dir_name[] = {"1550"}; 		// ��U�֐��̔g�� [nm]
//char *dir_name[] = {"1555", "1560", "1570", "1575"}; 		// ��U�֐��̔g�� [nm]
double lambda; 											// ��U�֐��̔g�� [m] (�v���O�������� dir_name ��double�����đ��)
double omega0; 											// ��U�֐��̊p���g�� [s^-1]

// �K�E�V�A���p���X
double sigma; 											// �L���蕝�����肷��萔
//static const double delta_omega = 0.1; 						// ���S���g���ŋK�i�������l�S��
//static const int Npeak = 3000; 								// �s�[�N�X�e�b�v��


/*****************************************************************************/
// Model�|���o���p�̒萔
/*****************************************************************************/
#define CLAD 0									// �N���b�h
#define CORE 1									// �R�A
#define GaInAsP	2								// GaInAsP
#define AIR_GaInAsP	3							// CLAD/GaInAsP
#define EXITATION 30							// ��U�_
#define OBSERVATION 20							// �ϑ��_
#define CIRCLE_REF_INDEX	CLAD				//�֐�mcircle�ŏ������ސ���
#define CIRCLE_REF_INDEX2	CLAD		//�֐�mcircle2�ŏ������ސ���
#define CIRCLE_REF_INDEX3	2				//�ꏊ�m�F�p�h�b�g�̐F�w��


/*****************************************************************************/
// �t�H�g�j�b�N�������f���p�錾
/*****************************************************************************/
//���������p�����[�^
struct PNUM {int X; int Y; }; 					//�~���̒��S���W��^����ϐ�.sankaku�Ŏg�p


/*****************************************************************************/
// �O���[�o���ϐ� (MPI�̒ʐM�ł̓O���[�o���ϐ����g�p���Ȃ��ƃG���[��������͗l)
/*****************************************************************************/

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION

// �d���E (XMAX+1�͕���v�Z�ł�"�̂肵��"���K�v�Ȃ��߁CYMAX+1, ZMAX+1�͑Ώ̋��E�����ŕK�v�Ȃ���)
double Ex[XMAX+1][YMAX+1][ZMAX+1];
double Ey[XMAX+1][YMAX+1][ZMAX+1];
double Ez[XMAX+1][YMAX+1][ZMAX+1];
double Hx[XMAX+1][YMAX+1][ZMAX+1];
double Hy[XMAX+1][YMAX+1][ZMAX+1];
double Hz[XMAX+1][YMAX+1][ZMAX+1];

// �U�d�� (YMAX+1, ZMAX+1�͑Ώ̋��E�����ŕK�v�Ȃ���)
double ALL_epsilonx[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];
double ALL_epsilony[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];
double ALL_epsilonz[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];
double epsilonx[XMAX][YMAX+1][ZMAX+1];
double epsilony[XMAX][YMAX+1][ZMAX+1];
double epsilonz[XMAX][YMAX+1][ZMAX+1];
int ALL_cell[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];		// ���ۂ�ALL_cell[XMAX_ALL+1][YMAX_ALL][ZMAX_ALL]�ŗǂ����C�v���O�����̓s����U�d�̂̔z��Ɠ��������ɂ���
int cell[XMAX][YMAX+1][ZMAX+1];							// ���ۂ�cell[XMAX][YMAX][ZMAX]�ŗǂ����C�v���O�����̓s����U�d�̂̔z��Ɠ��������ɂ���
double epsilon_xy[XMAX+1][YMAX+1], epsilon_yz[YMAX+1][ZMAX+1], epsilon_zx[XMAX+1][ZMAX+1];
double epsilon_zx2[XMAX+1][ZMAX+1];
int cell_xy[XMAX][YMAX], cell_yz[YMAX][ZMAX];

/*�z�����E�����K�p�̂Ƃ��Ɏg���d�E�̔z��*/
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

// ���E���z�o�͗p
double field_xy[XMAX][YMAX]; 	// Hz-field �̃t�@�C���o�� (�ʐ��������̎��E����)
double field_xz[XMAX][ZMAX]; 	// Hz-field �̃t�@�C���o�� (�ʐ��������̎��E����)
double field_yz[YMAX][ZMAX];	// Hz-field �̃t�@�C���o�� (�ʕ����̎��E����)

int x, y, z; 				// ���U���W
int xmax, ymax, zmax; 			// ��ԕ����̍ő�l
int xmax_all, ymax_all, zmax_all; //�����O�̍ő�l
//int n, Nmax; 				//n: ���ԃX�e�b�v�CNmax: ���ԃX�e�b�v�̍ő�l
int n; 					//���ԃX�e�b�v��
//int icut, jcut, kcut; 		//icut, jcut, kcut: �t�B�[���h�o�͂̂Ƃ��Ɏg��
//int PrintStat;
//int PrintEnd;
//int source_i1, source_i2;
//int source_j1, source_j2;
//int source_k;
int x_cen, y_cen, z_cen; 			//x_cen, y_cen, z_cen: ��͋�Ԃ̒��S���U���W
int x_model_cen, y_model_cen; 	//x_model_cen, y_model_cen:���f���̒��S���U���W
//char fname[30];
int irank, isize;

// �d���E�v�Z�Ɏg���萔�̐錾
static const double cnstHxyz = dt / MU0; 					//���E�̌v�Z�Ɏg���萔		cnstHxyz: Magnetic field calculation


///*��U�p�萔�̐錾*/
//double omega0; 	//omega0: ��U�֐��̊p���g��
//double sigma; 	//�K�E�V�A���p���X�̂��߂̍L���蕝�����肷��萔
//int Npeak; 		//�K�E�V�A���p���X�̃s�[�N�X�e�b�v��

/*���f���ݒ�Ɏg���萔�̐錾 -- �֐�modeling()�ȊO�ł��g���萔������̂ł����Ő錾����D*/
double b; //�f�B�X�N����
int bk; //�f�B�X�N�����̗��U�l

double w; 		//�x�����̔���
int wij; 		//�x�����̔����̗��U�l
double h; 		//�x���̍���
int hk; 			//�x���̍����̗��U�l
//double np; 		//�x���̋��ܗ�

//���o�̓p���[�̍ő�l�ƍŏ��l���L�^����ϐ�
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
/*�ϑ��_*/
int inputi, inputj, inputk;
int outputi, outputj, outputk;
double pyn_in, pyn_out;


/*�t�@�C���|�C���^�錾*/
//���f��			//Ez				//Hz
FILE *model_xy; 		FILE *fpfez_xy; 		FILE *fpfhz_xy;
FILE *model_yz; 		FILE *fpfez_yz; 		FILE *fpfhz_yz;
FILE *model_xz; 		FILE *fpfez_zx; 		FILE *fpfhz_zx;
FILE *model_xy2; 	FILE *fpfez2_xy; 	FILE *fpfhz2_xy;
FILE *allmodel_xy;
FILE *allmodel_yz1, *allmodel_yz4, *allmodel_yz7;

FILE *model_xy_1, *model_xy_2;

//�U�d��
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

FILE *fpparameter; 	//�v�Z�p�����[�^�ۑ��t�@�C��


#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
