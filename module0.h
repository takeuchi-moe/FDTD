/*****************************************************************************/
// ���̌v�Z�@�ł̉��͗̈�
/*****************************************************************************/
#define XMAX (INT_DIV(XMAX_ALL, NODE) + 1)
#define YMAX YMAX_ALL
#define ZMAX ZMAX_ALL


/*****************************************************************************/
// ���U�֐�
/*****************************************************************************/
double lambda; 											// ���U�֐��̔g�� [m] (�v���O�������� dir_name ��double�����đ���)
double omega0; 											// ���U�֐��̊p���g�� [s^-1]

// �K�E�V�A���p���X
double sigma; 											// �L���蕝�����肷���萔

/*****************************************************************************/
// Model�|���o���p�̒萔
/*****************************************************************************/
#define CLAD 0									// �N���b�h
#define CORE 1									// �R�A
#define GaInAsP	2								// GaInAsP
#define AIR_GaInAsP	3							// CLAD/GaInAsP
#define EXITATION 30							// ���U�_
#define OBSERVATION 20							// �ϑ��_
#define CIRCLE_REF_INDEX	CLAD				//�֐�mcircle�ŏ������ސ���
#define CIRCLE_REF_INDEX2	CLAD		//�֐�mcircle2�ŏ������ސ���
#define CIRCLE_REF_INDEX3	2				//�ꏊ�m�F�p�h�b�g�̐F�w��


/*****************************************************************************/
// �t�H�g�j�b�N�������f���p�錾
/*****************************************************************************/
//���������p�����[�^
struct PNUM {int X; int Y; }; 					//�~���̒��S���W���^�����ϐ�.sankaku�Ŏg�p


/*****************************************************************************/
// �O���[�o���ϐ� (MPI�̒ʐM�ł̓O���[�o���ϐ����g�p���Ȃ��ƃG���[���������͗l)
/*****************************************************************************/

#if _CALCULATION_TYPE == _PROPAGATION_CALCULATION

// �d���E (XMAX+1�͕����v�Z�ł�"�̂肵��"���K�v�Ȃ��߁CYMAX+1, ZMAX+1�͑Ώ̋��E�����ŕK�v�Ȃ���)
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
int ALL_cell[XMAX_ALL+1][YMAX_ALL+1][ZMAX_ALL+1];		// ���ۂ�ALL_cell[XMAX_ALL+1][YMAX_ALL][ZMAX_ALL]�ŗǂ����C�v���O�����̓s�����U�d�̂̔z���Ɠ��������ɂ���
int cell[XMAX][YMAX+1][ZMAX+1];							// ���ۂ�cell[XMAX][YMAX][ZMAX]�ŗǂ����C�v���O�����̓s�����U�d�̂̔z���Ɠ��������ɂ���
double epsilon_xy[XMAX+1][YMAX+1], epsilon_yz[YMAX+1][ZMAX+1], epsilon_zx[XMAX+1][ZMAX+1];
double epsilon_zx2[XMAX+1][ZMAX+1];
int cell_xy[XMAX][YMAX], cell_yz[YMAX][ZMAX], cell_xz[XMAX][ZMAX]; 


/*�z�����E�����K�p�̂Ƃ��Ɏg���d�E�̔z��*/
double Exn2y00[XMAX+1][ZMAX+1], 	Exn1y00[XMAX+1][ZMAX+1+1], 	Exn2y01[XMAX+1][ZMAX+1], 	Exn1y01[XMAX+1][ZMAX+1+1];
double Exn2z00[XMAX+1][YMAX+1], 	Exn1z00[XMAX+1][YMAX+1+1], 	Exn2z01[XMAX+1][YMAX+1], 	Exn1z01[XMAX+1][YMAX+1+1];


double Eyn2z00[XMAX+1][YMAX], 	Eyn1z00[XMAX+1][YMAX+1], 	Eyn2z01[XMAX+1][YMAX], 	Eyn1z01[XMAX+1][YMAX+1];
double Eyn2x00[YMAX][ZMAX+1], 	Eyn1x00[YMAX+1][ZMAX+1+1], 	Eyn2x01[YMAX][ZMAX+1], 	Eyn1x01[YMAX+1][ZMAX+1+1];
double Eyn2xm1[YMAX][ZMAX+1], 	Eyn1xm1[YMAX+1][ZMAX+1+1], 	Eyn2xm0[YMAX][ZMAX+1], 	Eyn1xm0[YMAX+1][ZMAX+1+1];

double Ezn2y00[XMAX+1][ZMAX], 	Ezn1y00[XMAX+1][ZMAX+1], 	Ezn2y01[XMAX+1][ZMAX], 	Ezn1y01[XMAX+1][ZMAX+1];
double Ezn2x00[YMAX+1][ZMAX], 	Ezn1x00[YMAX+1+1][ZMAX+1], 	Ezn2x01[YMAX+1][ZMAX], 	Ezn1x01[YMAX+1+1][ZMAX+1];
double Ezn2xm1[YMAX+1][ZMAX], 	Ezn1xm1[YMAX+1+1][ZMAX+1], 	Ezn2xm0[YMAX+1][ZMAX], 	Ezn1xm0[YMAX+1+1][ZMAX+1];

// ���E���z�o�͗p
double field_xy[XMAX][YMAX]; 	// Hz-field �̃t�@�C���o�� (�ʐ��������̎��E����)
double field_yz[YMAX][ZMAX];	// Hz-field �̃t�@�C���o�� (�ʕ����̎��E����)

int x, y, z; 				// ���U���W
int xmax, ymax, zmax; 			// ���ԕ����̍ő��l
int xmax_all, ymax_all, zmax_all; //�����O�̍ő��l
int n; 					//���ԃX�e�b�v��
int x_cen, y_cen, z_cen; 			//x_cen, y_cen, z_cen: ���͋��Ԃ̒��S���U���W
int x_model_cen, y_model_cen; 	//x_model_cen, y_model_cen:���f���̒��S���U���W
int irank, isize;

// �d���E�v�Z�Ɏg���萔�̐錾
static const double cnstHxyz = dt / MU0; 					//���E�̌v�Z�Ɏg���萔		cnstHxyz: Magnetic field calculation


///*���U�p�萔�̐錾*/

/*���f���ݒ��Ɏg���萔�̐錾 -- �֐�modeling()�ȊO�ł��g���萔�������̂ł����Ő錾�����D*/
double b; //�f�B�X�N����
int bk; //�f�B�X�N�����̗��U�l

double w; 		//�x�����̔���
int wij; 		//�x�����̔����̗��U�l
double h; 		//�x���̍���
int hk; 			//�x���̍����̗��U�l


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
FILE *fpfez2_xy; 	FILE *fpfhz2_xy;
FILE *allmodel_xy;


//�U�d��
FILE *fpepsilonx;
FILE *fpallepsilonx;
FILE *fpepsilony;
FILE *fpepsilonz;
//FILE *fpepsilony2;
//FILE *fpepsilonz2;
FILE *fpAllEpsilon;
FILE *fpEpsilon;


FILE *fpex; 		FILE *fphx;
FILE *fpey; 		FILE *fphy;
FILE *fpez; 		FILE *fphz;


FILE *fpparameter; 	//�v�Z�p�����[�^�ۑ��t�@�C��



#elif _CALCULATION_TYPE == _BAND_CALCULATION

#endif
