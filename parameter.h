/*****************************************************************************/
// �����v�Z���̌v�Z�@�̑䐔
/*****************************************************************************/
#if _FDTD
#define NODE 20
#else
#define NODE 1
#endif

// PCW�̍\��
#define PCW_Air_Or_SiO2	0		// 0:SiO2, 1:Air-brige
#define PCW_S1S3_Shift	0		// 0:3���ڊi�q�V�t�g�\��, 1:1,3���ڊi�q�V�t�g�\��

// �����v�Z���̌v�Z�@�̃����N
#define IRANK_MAX (NODE-1)		// �ő��l
#define IRANK_MIN 0				// �ŏ��l(���̐����̊��蓖�Ă��ꂽ�v�Z�@�Ɍ��ʂ�����)
#define ISIZE NODE				// �T�C�Y

#define TRUE 1
#define FALSE 0
/*****************************************************************************/
// �v�Z�p�����[�^		�P�ʂ�nm
// (�z���̒��`���v���v���Z�b�T�Ŏg�p���邽�߁C�}�N����)
// �e�l�̓Z���T�C�Y�̐����{�ɂȂ��Ă��邱��!!
/*****************************************************************************/

#if _PROGRAM_TEST

#if PCW_Air_Or_SiO2
/****************************** �{�ԗp(Air-brige) ******************************/

/*-------------------- CELL_SIZE:15nm --------------------*/
#define CELL_SIZE 15			// �Z���T�C�Y
#define PITCH 450				// PC �i�q�萔
#define PITCH_SHIFT_MAX 450 //480		// �i�q�萔�ω�PCW��PC�i�q�萔�̍ő��l

#define SLAB_HEIGHT 210		// �X���u��
#define CLAD_HEIGHT1 525	// �㕔�N���b�h����
#define CLAD_HEIGHT2 0		// �����N���b�h����
#define AIR_HEIGHT 0		// ���C�w����

#define RADIUS 155			// PC�̕W���~�E���a
#define SX1 0				// �`��(X)������1���ڊi�q�V�t�g��
#define SX3 90				// �`��(X)������3���ڊi�q�V�t�g��
#define SY 0				// ��(Y)�����̓��g�H�S�̊i�q�V�t�g��

#define EXCT_LEN 840					// ���U�_ (���f���̍��[�����̋���)
#define EXCT_OBSE_LEN 975				// ���U�_�����ϑ��ʂ̒��S�܂ł̋���
#define OBSE_WIRE_LEN 2540				// �ϑ��ʂ̒��S�����א����g�H�[�܂ł̋���
//#define OBSE_INTER (2 * PITCH)		// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define OBSE_INTER (PITCH)				// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 105)	// �o�ˍא����g�H�̒���(���ł��̐����̓Z���T�C�Y��NODE���̐����{�ɂ��邽�߂Ɏg�p)
#define WIRE_OUTPUT_OFFSET 0			// �o�ˍא����g�H�̃X���u�I�[�̒���
#define WIRE_WID_OFFSET 0				// �א����̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_SiSLAB_TERMINATION_LEN 255	// PCW����CORE�X���u�I�[�̒���
#define PCW_SiSLAB_OFFSET 0				// PCW�c��CORE�X���u�̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_WIDTH_CHIRP 0				// PCW���̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)

// �����ȉ��͎�����
#define NORM_PCW_PER 0				// �ʏ�PCW������
#define CHIRP_3RD_LS_PER 0			// 3���ڊi�q�V�t�g�ʃ`���[�vLSPCW������
#define PITCH_SHIFT_PER 5			// �i�q�萔�ω�PCW�̎�����
#define LSPCW_SHIFT_DESCRETE FALSE			// �i�q�萔�ω�PCW�̂Ƃ��C�͂��߂���LSPCW�ɂ����ꍇ��FALSE�@���Ȃ��ꍇ��TRUE
#define PITCH_SHIFT_CHIRP_PER 0		// �i�q�萔�ω��`���[�vPCW�̎�����
#define LSPCW_PER 20				// LSPCW������
#define PCW_WID 6					// PCW�̗���
// �����ȏ��͎�����
#define LSPCW_ROW 3					// ���g�H���琔���ĉ����ڂ̊i�q�_���`�������ɃV�t�g�����邩�D

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// �ʏ�PCW��
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW��
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// �o�ˍא����g�H�̃X���u�I�[�̒���

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// ���ˊϑ��ʂ̒��S���W (���f���̍��[�����̋���)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// �o�ˊϑ��ʂ̒��S���W (���f���̉E�[�����̋���)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// ���ˍא���
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// �o�ˍא���
#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW��
/*-------------------- CELL_SIZE:15nm --------------------*/

/*-------------------- CELL_SIZE:21nm --------------------*/
//#define CELL_SIZE 21			// �Z���T�C�Y
//#define PITCH 441				// PC �i�q�萔
//#define PITCH_SHIFT_MAX 441 //480		// �i�q�萔�ω�PCW��PC�i�q�萔�̍ő��l
//
//#define SLAB_HEIGHT 210		// �X���u��
//#define CLAD_HEIGHT1 525	// �㕔�N���b�h����
//#define CLAD_HEIGHT2 0		// �����N���b�h����
//#define AIR_HEIGHT 0		// ���C�w����
//
//#define RADIUS 147			// PC�̕W���~�E���a
//#define SX3 105				// �`��(X)������3���ڊi�q�V�t�g��
//#define SY 0				// ��(Y)�����̓��g�H�S�̊i�q�V�t�g��
//
//#define EXCT_LEN 840					// ���U�_ (���f���̍��[�����̋���)
//#define EXCT_OBSE_LEN 966				// ���U�_�����ϑ��ʂ̒��S�܂ł̋���
//#define OBSE_WIRE_LEN 2541				// �ϑ��ʂ̒��S�����א����g�H�[�܂ł̋���
////#define OBSE_INTER (2 * PITCH)		// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
//#define OBSE_INTER (PITCH)				// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
//#define WIRE_OUTPUT_LEN (EXCT_LEN + 0)	// �o�ˍא����g�H�̒���(���ł��̐����̓Z���T�C�Y��NODE���̐����{�ɂ��邽�߂Ɏg�p)
//#define WIRE_OUTPUT_OFFSET 0			// �o�ˍא����g�H�̃X���u�I�[�̒���
//#define WIRE_WID_OFFSET 0				// �א����̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
//#define PCW_SiSLAB_TERMINATION_LEN 252	// PCW����CORE�X���u�I�[�̒���
//#define PCW_SiSLAB_OFFSET 0				// PCW�c��CORE�X���u�̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
//#define PCW_WIDTH_CHIRP 0				// PCW���̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
//
//// �����ȉ��͎�����
//#define NORM_PCW_PER 0				// �ʏ�PCW������
//#define CHIRP_3RD_LS_PER 0			// 3���ڊi�q�V�t�g�ʃ`���[�vLSPCW������
//#define PITCH_SHIFT_PER 5			// �i�q�萔�ω�PCW�̎�����
//#define PITCH_SHIFT_CHIRP_PER 0		// �i�q�萔�ω��`���[�vPCW�̎�����
//#define LSPCW_PER 20				// LSPCW������
//#define PCW_WID 6					// PCW�̗���
//// �����ȏ��͎�����
//
//#define LSPCW_ROW 3					// ���g�H���琔���ĉ����ڂ̊i�q�_���`�������ɃV�t�g�����邩�D
//
//#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// �ʏ�PCW��
//#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
//#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW��
//#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// �o�ˍא����g�H�̃X���u�I�[�̒���
//
//#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// ���ˊϑ��ʂ̒��S���W (���f���̍��[�����̋���)
//#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// �o�ˊϑ��ʂ̒��S���W (���f���̉E�[�����̋���)
//
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// ���ˍא���
//#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// �o�ˍא���
//#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW��
/*-------------------- CELL_SIZE:21nm --------------------*/

/****************************** �{�ԗp(Air-brige) ******************************/
#else

/****************************** �{�ԗp(SiO2) ******************************/

#if PCW_S1S3_Shift

/****************************** 1,3���ڊi�q�V�t�g�\�� ******************************/
#define CELL_SIZE 15			// �Z���T�C�Y
#define PITCH 405				// PC �i�q�萔
#define PITCH_SHIFT_MAX 495 //480		// �i�q�萔�ω�PCW��PC�i�q�萔�̍ő��l

#define SLAB_HEIGHT 210		// �X���u��
#define CLAD_HEIGHT1 750	// �㕔�N���b�h����
//#define CLAD_HEIGHT1 990	// �㕔�N���b�h����
#define CLAD_HEIGHT2 0		// �����N���b�h����
#define AIR_HEIGHT 0		// ���C�w����

#define RADIUS 110			// PC�̕W���~�E���a
//#define SX3 90				// �`��(X)������3���ڊi�q�V�t�g��
//#define SX1 0				// �`��(X)������3���ڊi�q�V�t�g��
//#define SX3 90				// �`��(X)������3���ڊi�q�V�t�g��
//#define SX1 120				// �`��(X)������1���ڊi�q�V�t�g��
#define SX3 0				// �`��(X)������3���ڊi�q�V�t�g��
#define SX1 0				// �`��(X)������1���ڊi�q�V�t�g��
#define SY 15				// ��(Y)�����̓��g�H�S�̊i�q�V�t�g��(�Ώ̋��E���g�p���Ă����̂Ŏ��ۂɂ͂��̔{)

#define EXCT_LEN 840					// ���U�_ (���f���̍��[�����̋���)
#define EXCT_OBSE_LEN 975				// ���U�_�����ϑ��ʂ̒��S�܂ł̋���
#define OBSE_WIRE_LEN 2540				// �ϑ��ʂ̒��S�����א����g�H�[�܂ł̋���
//#define OBSE_INTER (2 * PITCH)		// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define OBSE_INTER (PITCH)				// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 60)	// �o�ˍא����g�H�̒���(���ł��̐����̓Z���T�C�Y��NODE���̐����{�ɂ��邽�߂Ɏg�p)
#define WIRE_OUTPUT_OFFSET 0			// �o�ˍא����g�H�̃X���u�I�[�̒���
#define WIRE_WID_OFFSET 0				// �א����̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_SiSLAB_TERMINATION_LEN 240	// PCW����CORE�X���u�I�[�̒���
#define PCW_SiSLAB_OFFSET 0				// PCW�c��CORE�X���u�̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_WIDTH_CHIRP 0				// PCW���̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)

// �����ȉ��͎�����
#define NORM_PCW_PER 0				// �ʏ�PCW������
#define CHIRP_3RD_LS_PER 0			// 3���ڊi�q�V�t�g�ʃ`���[�vLSPCW������
#define PITCH_SHIFT_PER 0			// �i�q�萔�ω�PCW�̎�����
#define PITCH_SHIFT_CHIRP_PER 5	// �i�q�萔�ω��`���[�vPCW�̎�����
#define LSPCW_PER 5				// LSPCW������
#define PCW_WID 6					// PCW�̗���
// �����ȏ��͎�����

//#define LSPCW_ROW 3					// ���g�H���琔���ĉ����ڂ̊i�q�_���`�������ɃV�t�g�����邩�D

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// �ʏ�PCW��
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW��
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// �o�ˍא����g�H�̃X���u�I�[�̒���

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// ���ˊϑ��ʂ̒��S���W (���f���̍��[�����̋���)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)	// �o�ˊϑ��ʂ̒��S���W (���f���̉E�[�����̋���)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)			// ���ˍא���
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)		// �o�ˍא���
//#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW��

/****************************** 1,3���ڊi�q�V�t�g�\�� ******************************/

#else

/****************************** 3���ڊi�q�V�t�g�\�� ******************************/
#define CELL_SIZE 21			// �Z���T�C�Y   15
#define PITCH 399				// PC �i�q�萔   405
#define PITCH_SHIFT_MAX 399 //480		// �i�q�萔�ω�PCW��PC�i�q�萔�̍ő��l
#define PITCH_SHIFT_MAX2 399 //���`���[�v�̍ŏ��l �����݂̂�����

#define SLAB_HEIGHT 210		// �X���u��
#define CLAD_HEIGHT1 1995	// �㕔�N���b�h����   750
//#define CLAD_HEIGHT1 990	// �㕔�N���b�h����
#define CLAD_HEIGHT2 0		// �����N���b�h����
#define AIR_HEIGHT 0		// ���C�w����

#define RADIUS 110			// PC�̕W���~�E���a      120
#define SX3 84				// �`��(X)������3���ڊi�q�V�t�g��(SX2,SX4=0�łȂ��Ǝg���Ȃ��v���P!!main��1380�s��)
#define SX1 0				// �`��(X)������1���ڊi�q�V�t�g��
#define SX2 0				// �`��(X)������2���ڊi�q�V�t�g��(SX3,SX4=0�łȂ��Ǝg���Ȃ��v���P!!main��1380�s��)
#define SX4 0				//�`��(X)������4���ڊi�q�V�t�g��(SX2,SX3=0�łȂ��Ǝg���Ȃ��v���P!!main��1380�s��)
#define SY -588				//-588 //������(Y)�����̓��g�H�S�̊i�q�V�t�g��(�Ώ̋��E���g�p���Ă����̂Ŏ��ۂɂ͂��̔{)�@-588����+�ŕ��L

#define EXCT_LEN 840			//840		// �����U�_ (���f���̍��[�����̋���) �����͊��{�ς��Ȃ��D
#define EXCT_OBSE_LEN 3000-420		//975 		// �������U�_�����ϑ��ʂ̒��S�܂ł̋���  ���������{�ς��Ȃ�  3000//15/12/25
//OBSE_WIRE_LEN��+a�����Ƃ��ɂ͂�����-a����
#define OBSE_WIRE_LEN 5815		//2540(121 �Z��) 2750(131�Z��) 5080(242�Z��)//  5500(262�Z��)2/23 //5815(297�Z��) 35�Z�����₵�� 2/24  �����ɒ��ډ�(by�n���搶) ��(�o��)�ϑ��ʂ̒��S�����א����g�H�[�܂ł̋���
//������+a�����ƁC���j�^��PCW�[�̋�����+a�C���j�^�Əo�˒[�̋�����-2a�������D�o�˂̍א�������-a�ɂȂ�.���˃��j�^��PCW�[�Ƃ̋�����+a
//#define OBSE_INTER (2 * PITCH)		// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define OBSE_INTER (PITCH)				// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define WIRE_OUTPUT_LEN (EXCT_LEN + 45)	// �o�ˍא����g�H�̒���(���ł��̐����̓Z���T�C�Y��NODE���̐����{�ɂ��邽�߂Ɏg�p)
#define WIRE_OUTPUT_OFFSET 0			// �o�ˍא����g�H�̃X���u�I�[�̒���
#define WIRE_WID_OFFSET 168 //168				// �א����̃I�t�Z�b�g�ʂ̔���(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define WIRE_WID_OFFSET_OUT 168//168
#define PCW_SiSLAB_TERMINATION_LEN 255	//255 // ��PCW����CORE�X���u�I�[�̒����@+�ɂ����Ɖ~�E�̂ݒ��S�����ɂ����� (�����̒P�ʂ�nm)
#define PCW_SiSLAB_OFFSET 0 //0				// ��PCW�c��CORE�X���u�̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)2
//�`�������V�t�g
#define PCW_WIDTH_CHIRP 168//168				// PCW���̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_WIDTH_CHIRP_OUT 168//168
#define PCW_WIDTH_Para 1

//�����g�H���`���[�v�ʂ̃p�����[�^�S��(4��)��0�ɂ����ƃo�O���������̂ŕ֋X�I�ɑS��1�ɂ���

// �����ȉ��͎�����
#define NORM_PCW_PER 0				// 160722 ������0�Œ� �ʏ�PCW������
#define CHIRP_3RD_LS_PER 0			// 160722 ������0�Œ� 3���ڊi�q�V�t�g�ʃ`���[�vLSPCW�������D 2~4���ڂ܂őΉ��D���V�t�g��/������<cellsize���ƃ`���[�v���Ȃ�
#define CHIRP_2ND_LS_PER 7//5			//�����g�H���`���[�v�ƃV�t�g�ʃ`���[�v�𓯎��ɍs���ۂ̃`���[�vLSPCW�������DPITCH_SHIFT_PER���菬�����Ȃ��ƃ_��?�@2,3���ڂ̂ݑΉ��D���V�t�g��/������<cellsize���ƃ`���[�v���Ȃ�
#define CHIRP_2ND_LS_PER_OUT 7//5 ��
#define PITCH_SHIFT_PER 7//5				//�������g�H���`���[�v�̎������@16/7/22 ���ꂪ7�ȏゾ�ƃ_���� (6��OK)
#define PITCH_SHIFT_PER_OUT 7//5 ����
#define PITCH_SHIFT_CHIRP_PER2 7//5		// ���i�q�萔�ω��`���[�vPCW�̎�����(���g�H���`���[�v�Ɠ���)
#define PITCH_SHIFT_CHIRP_PER2_OUT 7//��5
#define LSPCW_SHIFT_DESCRETE FALSE			// �i�q�萔�ω�PCW�̂Ƃ��C�͂��߂���LSPCW�ɂ����ꍇ��FALSE�@���Ȃ��ꍇ��TRUE
#define PITCH_SHIFT_CHIRP_PER 0		// �i�q�萔�ω��`���[�vPCW�̎�����
#define LSPCW_PER 16 //LSPCW������
#define PCW_WID 8 //8 PCW�̗��� //SIWG�ł�1
// �����ȏ��͎�����
				//
//#define LSPCW_ROW 3					// ���g�H���琔���ĉ����ڂ̊i�q�_���`�������ɃV�t�g�����邩�D

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)				// �ʏ�PCW��
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
#define LSPCW_LEN (LSPCW_PER * PITCH)					// LSPCW��
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)	// �o�ˍא����g�H�̃X���u�I�[�̒���

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)			// ���ˊϑ��ʂ̒��S���W (���f���̍��[�����̋���)
//#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN - 210)	// ���@�o�ˊϑ��ʂ̒��S���W (���f���̉E�[�����̋���)�@�����͎��ۂ̒���(nm)
//���̈ʒu�͂����炭�}�ɂ͕\�������Ȃ��@?intObseLenPart4�Ƃ͌������t
//PS.15/12/22 �����͋@�\���ĂȂ�



#define WIRE_LEN1 (OBSE_LEN1 + 3696)			// �����ˍא����@//2541 ���� //2961�@�����ɒ���2/23 //�����ɒ���2 3696 2/24
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)�@//���̃v���O����
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)

#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN) // ���� �o�ˍא��� ������nm�@�o�˒[�����o�˃��j�^�̋���
//���j�^�[�������ɂ����ꍇ�̓I�t�Z�b�g137*21
//#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW��

/****************************** 3���ڊi�q�V�t�g�\�� ******************************/

#endif
/****************************** �{�ԗp(SiO2) ******************************/
#endif

#else

/****************************** ���͗̈挈���p ******************************/

#define CELL_SIZE 35		// �Z���T�C�Y
#define PITCH 455			// PC �i�q�萔
#define PITCH_SHIFT_MAX 455 // �i�q�萔�ω�PCW��PC�i�q�萔�̍ő��l

#define SLAB_HEIGHT 210		// �X���u��
#define CLAD_HEIGHT1 525	// �㕔�N���b�h����
#define CLAD_HEIGHT2 0		// �����N���b�h����
#define AIR_HEIGHT 0		// ���C�w����

#define RADIUS 160			// PC�̕W���~�E���a
#define SX3 105				// �`��(X)������3���ڊi�q�V�t�g��
#define SX1 0				// �`��(X)������3���ڊi�q�V�t�g��
#define SY 0				// ��(Y)�����̓��g�H�S�̊i�q�V�t�g��

#define EXCT_LEN 840					// ���U�_ (���f���̍��[�����̋���)
#define EXCT_OBSE_LEN 980				// ���U�_�����ϑ��ʂ̒��S�܂ł̋���
#define OBSE_WIRE_LEN 2520				// �ϑ��ʂ̒��S�����א����g�H�[�܂ł̋���
#define OBSE_INTER (2 * PITCH)			// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
#define WIRE_OUTPUT_LEN EXCT_LEN		// �o�ˍא����g�H�̒���
#define WIRE_OUTPUT_OFFSET 0			// �o�ˍא����g�H�̃X���u�I�[�̒���
#define WIRE_WID_OFFSET 0				// �א����̃I�t�Z�b�g��(�v���X�������L)
#define PCW_SiSLAB_TERMINATION_LEN 105	// PCW����CORE�X���u�I�[�̒���
#define PCW_SiSLAB_OFFSET 0				// PCW�c��CORE�X���u�̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)
#define PCW_WIDTH_CHIRP 0				// PCW���̃I�t�Z�b�g��(�v���X�������L) (�����̏ꍇ�ɂ́C1�Z�����傫���ݒ肵�Ȃ��Ɗۂߌ덷����)

// �����ȉ��͎�����
#define NORM_PCW_PER 0				// �ʏ�PCW������
#define CHIRP_3RD_LS_PER 0			// 3���ڊi�q�V�t�g�ʃ`���[�vLSPCW������
#define PITCH_SHIFT_PER 5			// �i�q�萔�ω�PCW�̎�����
#define PITCH_SHIFT_CHIRP_PER 0		// �i�q�萔�ω��`���[�vPCW�̎�����
#define LSPCW_PER 20				// LSPCW������
#define PCW_WID 8					// PCW�̗���
// �����ȏ��͎�����

#define LSPCW_ROW 3					// ���g�H���琔���ĉ����ڂ̊i�q�_���`�������ɃV�t�g�����邩�D

#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)		// �ʏ�PCW��
#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
#define LSPCW_LEN (LSPCW_PER * PITCH)			// LSPCW��
#define WIRE_OUTPUT_OFFSET_PER INT_DIV(WIRE_OUTPUT_OFFSET, CELL_SIZE)			// �o�ˍא����g�H�̃X���u�I�[�̒���
// �����ȏ��͎�����

#define OBSE_LEN1 (EXCT_LEN + EXCT_OBSE_LEN)	// ���ˊϑ��ʂ̒��S���W (���f���̍��[�����̋���)
#define OBSE_LEN5 (WIRE_OUTPUT_OFFSET + EXCT_OBSE_LEN)		// �o�ˊϑ��ʂ̒��S���W (���f���̉E�[�����̋���)

#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)						// ���ˍא���
#define WIRE_LEN2 (OBSE_WIRE_LEN + WIRE_OUTPUT_LEN)					// �o�ˍא���
#define PCW_LEN (NORM_PCW_LEN * 2 + CHIRP_3RD_LS_LEN * 2 + LSPCW_LEN)	// PCW��

/****************************** ���͗̈挈���p ******************************/


/****************************** �����m�F�p ******************************/

//#define CELL_SIZE 100		// �Z���T�C�Y
//#define PITCH 400			// PC �i�q�萔
//
//#define SLAB_HEIGHT 200		// �X���u��
//#define CLAD_HEIGHT1 700	// �㕔�N���b�h����
//#define CLAD_HEIGHT2 0		// �����N���b�h����
//#define AIR_HEIGHT 0		// ���C�w����
//
//#define RADIUS 100			// PC�̕W���~�E���a
//#define SX3 100				// �`��(X)������3���ڊi�q�V�t�g��
//#define SY 0				// ��(Y)�����̓��g�H�S�̊i�q�V�t�g��
//
//#define NORM_PCW_PER 5			// �ʏ�PCW������
//#define CHIRP_3RD_LS_PER 1			// �`���[�vLSPCW������
//#define LSPCW_PER 5				// LSPCW������
//#define PCW_WID 6					// PCW�̗���
//#define LSPCW_ROW 3					// ���g�H���琔���ĉ����ڂ̊i�q�_���`�������ɃV�t�g�����邩�D
//#define NORM_PCW_LEN (NORM_PCW_PER * PITCH)		// �ʏ�PCW��
//#define CHIRP_3RD_LS_LEN (CHIRP_3RD_LS_PER * PITCH)		// �`���[�vLSPCW��
//#define LSPCW_LEN (LSPCW_PER * PITCH)			// LSPCW��
//
//#define EXCT_LEN 800					// ���U�_ (���f���̍��[�����̋���)
//#define EXCT_OBSE_LEN 1900				// ���U�_�����ϑ��ʂ̒��S�܂ł̋���
//#define OBSE_WIRE_LEN EXCT_OBSE_LEN		// �ϑ��ʂ̒��S�����א����g�H�܂ł̋���
//#define OBSE_INTER (2 * PITCH)			// �ϑ��ʂ̒���(���ߗ��C���˗������߂邽�߂Ɏg�p)
//#define OBSE_LEN1 (EXCT_LEN + OBSE_INTER / 2 + EXCT_OBSE_LEN)	// ���ˊϑ��ʂ̒��S���W (���f���̍��[�����̋���)
//#define OBSE_LEN5 (EXCT_LEN + OBSE_INTER / 2)		// �o�ˊϑ��ʂ̒��S���W (���f���̉E�[�����̋���)
//
//#define WIRE_LEN1 (OBSE_LEN1 + OBSE_WIRE_LEN)	// ���ˍא���
//#define WIRE_LEN2 (OBSE_LEN5 + OBSE_WIRE_LEN)	// �o�ˍא���
//#define WIRE_WID_OFFSET 0		// �א����̃I�t�Z�b�g��

/****************************** �����m�F�p ******************************/

#endif

/*****************************************************************************/
// �S���͗̈�
/*****************************************************************************/
#if _PROGRAM_TEST

#if PCW_Air_Or_SiO2
/*-------------------- CELL_SIZE:15nm --------------------*/
#define XMAX_ALL 1440	// SiO2 1340, Air 1440
#define YMAX_ALL 183	// SiO2 163 Air 183
#define ZMAX_ALL 42
/*-------------------- CELL_SIZE:15nm --------------------*/
#else
#if PCW_S1S3_Shift
/*-------------------- CELL_SIZE:15nm --------------------*/
#define XMAX_ALL 1725	// SiO2 1340
#define YMAX_ALL 180	// PCW_WID:6/163  PCW_WID:8/209  PCW_WID:10/255
#define ZMAX_ALL 57		// CLAD_HEIGHT1:525/42  CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73
//#define ZMAX_ALL 73	// CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73
/*-------------------- CELL_SIZE:15nm --------------------*/
#else
/*-------------------- CELL_SIZE:15nm --------------------*/
//��15nm�͋C�ɂ��Ȃ��ėǂ��D
#define XMAX_ALL 1335 //948(2015/11/11) //1250(2015/11/19) 1250-165=1085(2015/12/16) 1250-165-16=1069(2015/12/24) 1069+151-55=1165(2015/12/25)���˒��ډ�  1220 �Ƃ�����1206�ɂȂ�(2065/2/23)���˒��ډ�1
//1222 �� 1220 //1240 �� 1235 //1260 �� 1235 //1270 �� 1235 1240�����傫�����Ă��Ӗ��Ȃ��̂Ŋ��{1240
//1235+(6*400*2)/21=1235+230=1465
//1225 BOUNDARYLINE 41, 42(16/10/21)
//1245(1236) (16/10/24)
#define YMAX_ALL 172	//10���ł�145+32(���̏ꍇXMAX��50���x���炷�K�v����)�D8���ł�145.���ڂł�172
#define ZMAX_ALL 100	//15/12/16 70(�e�ʕs���̂���) �]��100
//#define ZMAX_ALL 73	// CLAD_HEIGHT1:750/57  CLAD_HEIGHT1:990/73
/*-------------------- CELL_SIZE:15nm --------------------*/
#endif
#endif

/*-------------------- CELL_SIZE:21nm --------------------*/
//#define XMAX_ALL 1010 // SiO2 1340, Air 695
//#define YMAX_ALL 127 //SiO2 163 Air 183
//#define ZMAX_ALL 30
/*-------------------- CELL_SIZE:21nm --------------------*/

#else

#define XMAX_ALL 616
#define YMAX_ALL 95
#define ZMAX_ALL 18

//#define XMAX_ALL 950
//#define YMAX_ALL 114
//#define ZMAX_ALL 38

//#define XMAX_ALL 162
//#define YMAX_ALL 19
//#define ZMAX_ALL 8

#endif

/*****************************************************************************/
// �Z���T�C�Y [m]
/*****************************************************************************/
static const double dblCellSize = CELL_SIZE * 1e-9;
static const double dx = dblCellSize;
static const double dy = dblCellSize;
static const double dz = dblCellSize;
static const double inv_dx = 1/dx;
static const double inv_dy = 1/dy;
static const double inv_dz = 1/dz;

/*****************************************************************************/
// ���ԃX�e�b�v (�N�[�����g�̈��������Ȃǂɒ���)
/*****************************************************************************/
#if _FDTD

#if _PROGRAM_TEST
#if PCW_Air_Or_SiO2
/*-------------------- CELL_SIZE:15nm --------------------*/
static const double dt = 28e-18; 			// ���ԃX�e�b�v[s]
static const int Nmax = 150000; 				// �ŏI���ԃX�e�b�v
/*-------------------- CELL_SIZE:15nm --------------------*/
#else
/*-------------------- CELL_SIZE:15nm --------------------*/ //��
static const double dt = 38e-18; 			// ���ԃX�e�b�v[s]
static const int Nmax = 150000; //150000 				// ���ŏI���ԃX�e�b�v
//static const int Nmaxp = 5000;  				// �����ŏI���ԃX�e�b�v
/*-------------------- CELL_SIZE:15nm --------------------*/
#endif

/*-------------------- CELL_SIZE:21nm --------------------*/
//static const double dt = 39e-18; 			// ���ԃX�e�b�v[s]
//static const int Nmax = 130000; 				// �ŏI���ԃX�e�b�v
//static const int Nmax = 250000; 				// �ŏI���ԃX�e�b�v
/*-------------------- CELL_SIZE:21nm --------------------*/


static const int Ncut = 50000; //50000				// �����ԃX�e�b�v�����\���������Ԋu

static const int Tcut = 1000; 	//1000			//�@�����p���[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v  (�ŏI�v�Z�X�e�b�v�����̍�)
//static const int Fcut = 500; 				// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)
static const int Fcut = 200; 				//200 �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)
//���ԃX�e�b�v����10000���ƁC
#else

static const double dt = 67e-18; 			// ���ԃX�e�b�v[s]

/******************** �������Z�� (<<1min) ********************/
//static const int Ncut = 5; 				// ���ԃX�e�b�v�����\���������Ԋu
//static const int Tcut = 30; 			// �G�l���M�[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v
//static const int Fcut = 99; 			// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)
//static const int Nmax = 100; 			// �ŏI���ԃX�e�b�v

/******************** �����Ȃ��Z�� (<5min) ********************/
static const int Ncut = 500; 			// ���ԃX�e�b�v�����\���������Ԋu
static const int Tcut = 200; 			// �G�l���M�[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v
static const int Fcut = 200; 			// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)
static const int Nmax = 1000; 			// �ŏI���ԃX�e�b�v

/******************** �����Z�� (<10min) ********************/
//static const int Ncut = 500; 				// ���ԃX�e�b�v�����\���������Ԋu
//static const int Tcut = 3000; 			// �G�l���M�[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v
//static const int Fcut = 500; 				// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)
//static const int Nmax = 10000; 			// �ŏI���ԃX�e�b�v

/******************** ����********************/
//static const int Ncut = 20000; 				// ���ԃX�e�b�v�����\���������Ԋu
//static const int Nmax = 100000; 			// �ŏI���ԃX�e�b�v
//static const int Tcut = Nmax - 100; 			// �G�l���M�[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v
//static const int Fcut = 200; 				// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)

#endif

#else
static const double dt = 2.8e-17; 			// ���ԃX�e�b�v[s]
static const int Ncut = 5; 				// ���ԃX�e�b�v�����\���������Ԋu
static const int Tcut = 30; //30			// ���G�l���M�[�̕��ς̎Z�o���J�n���鎞�ԃX�e�b�v

static const int Fcut = 30; 			// �t�B�[���h���o�͂��鎞�ԃX�e�b�v�� (�ŏI�v�Z�X�e�b�v�����̍�)
static const int Nmax = 1; 			// �ŏI���ԃX�e�b�v
#endif

//PoyntingPower

static const int Ncheck = 10; //10					// �������m�F�p�̃t�B�[���h���o�͂��鎞�ԃX�e�b�v
// �����͏��߂ɂ�������
static const int Ncutfield = Ncut; 			// �t�B�[���h���o�͂��鎞�ԃX�e�b�v
//static const int Ncutfield2 = 10; 			// �������Ԃł̃t�B�[���h���o�͂��鎞�ԃX�e�b�v�Ԋu
static const int Ncutfield2 = 5; 			// ���������Ԃł̃t�B�[���h���o�͂��鎞�ԃX�e�b�v�Ԋu
//������149800�`150000�܂�5���݂ŏo�͂���

/*****************************************************************************/
// ������[MKSA�n]
/*****************************************************************************/
#define PI 3.141592
#define C0 2.997924e8			// �^�󒆂̌��� [m/s]
#define EPSILON0 8.854e-12		// �^�󒆂̗U�d�� [F/m]
#define MU0 (PI*(4.0e-7))		// �^�󒆂̓����� [N/A^2]


/*****************************************************************************/
// �X���u�Ɖ��͗̈�
/*****************************************************************************/
static const double dblSlabHeig = SLAB_HEIGHT * 1.0e-9; 				// �X���u��
static const double dblCladHeight1 = CLAD_HEIGHT1 * 1.0e-9; 				// �㕔�N���b�h�w��
static const int	intSlabHeigPer = INT_DIV (dblSlabHeig/2.0, dblCellSize); 				//�X���u����(�Ώ̋��E�������g�p���Ă����̂ŁC�X���u����1/2)
static const int	intCladHeight1 = INT_DIV (CLAD_HEIGHT1, CELL_SIZE); 				//�㕔�N���b�h����
static const int	air_hc = 	(int)((0.0e-6*1e10)/(dz*1e10)); 					//�N���b�h�㕔���C�w��
static const int	intSlabCen = air_hc + intCladHeight1 + intSlabHeigPer; 			//�����w���S�Z��(+1��[intSlabHeigPer]����Z�����̂Ƃ��ɒ����l�ɂ��邽��)


/*****************************************************************************/
// �t�H�g�j�b�N�������g�H
/*****************************************************************************/
static const double dblPitchCellComp = dblCellSize / 2.0; 			// �i�q�萔�̊ۂ�
static const double dblPitch = PITCH * 1e-9; 					// �~�E�i�q�萔

static const double dblRadius = RADIUS * 1e-9; 					// �~�E���a
static const double dblRadius2 = 0.18e-6; 					// ���~�E���a
static const double dblRadius3 = 0.12e-6; 					// ���˕��~�E�̔��a
static const double dblRadius4 = 0.04e-6; 					// �ꏊ�m�F�p�h�b�g�̔��a
static const double dblRadius5 = 0.06e-6; 						//�����~�E�����[�~�E���a
static const double dblRadius6 = 0.10e-6; 						//�����~�E�����[�~�E���a
static const double dblRadius7 = 0.14e-6; 						//�����~�E�����[�~�E���a
static const double dblRadius8 = 0.16e-6; 						//�����~�E�����[�~�E���a

static const double dblRadius9 = 0.078e-6; 						//�e�[�p�p�~�E�T�C�Y 16/08/27
static const double dblRadius10 = 0.08e-6; 						//�e�[�p�p�~�E�T�C�Y 16/08/27
static const double dblRadius11 = 0.08e-6; 					//�e�[�p�p�~�E�T�C�Y 16/08/27
static const double dblRadius12 = 0.084e-6; 						//�e�[�p�p�~�E�T�C�Y 16/08/27
static const double dblRadius13 = 0.084e-6; 					//�e�[�p�p�~�E�T�C�Y 16/08/27
static const double dblRadius14 = 0.10e-6; 					//�e�[�p�p�~�E�T�C�Y 16/08/29 168
static const double dblRadius15 = 0.10e-6; 					//�e�[�p�p�~�E�T�C�Y 16/08/29 180

static const double dblDiamter = 2.0 * dblRadius; 				// �~�E���a

static const int intPitchX = INT_DIV(PITCH, CELL_SIZE); 		// �i�q�萔�̃Z���T�C�Y(X����)
static const int intPitchY = (INT_DIV((PITCH * sqrt(3.0) / 2 + 0.5), CELL_SIZE)); 		// �i�q�萔�̃Z���T�C�Y(Y����)	+0.5�͎l�̌ܓ��̂���
static const int intRadius = INT_DIV(RADIUS, CELL_SIZE); 		// �~�E���a


//static const int intPcwWid = 9; 		// �~�E�s��(�w������)(Ex.13)
//static const int intPcwLen = 51; 		// �~�E����(���g�H����)(Ex.10)
//static const int intPcwStartX = 8; 		// �~���z�u�̊J�nX���W(���g�H�����C�Z��������)(Ex.105)
//static const int intPcwStartY = 8; 		// �~���z�u�̊J�nY���W(�������C�Z��������)(Ex.4)
//static const int PCmargin = 6; 		//�t�H�g�j�b�N�����̈��}�[�W���D���̒��������ŊO�~�E�̒��S�����`�������ɃX�y�[�X���݂���(�Z�������́C"0" = �E�ʊԊu0.08um)

static const int intNormPcwPer = NORM_PCW_PER; 				// PCW����
static const int intPitchShiftPcwPer = PITCH_SHIFT_PER;				// ���i�q�萔�ω�PCW�̎�����
static const int intPitchShiftPcwPerOut = PITCH_SHIFT_PER_OUT; //��
static const int intPitchShiftChirpPcwPer = PITCH_SHIFT_CHIRP_PER;	// �i�q�萔�ω��`���[�vPCW�̎�����

static const int intChirp3rdLsPer = CHIRP_3RD_LS_PER; 		// �`���[�vLSPCW����
static const int intChirp2ndLsPer = CHIRP_2ND_LS_PER; 		// �`���[�vLSPCW����
static const int intChirp2ndLsPerOut = CHIRP_2ND_LS_PER_OUT;
static const int intLspcwPer = LSPCW_PER; 					// LSPCW����
static const int intPcwPer = 2 * (intNormPcwPer + intChirp3rdLsPer + intPitchShiftChirpPcwPer) + intPitchShiftPcwPer + intPitchShiftPcwPerOut + intLspcwPer; 	// ���SPCW����
static const int intSx4Per = INT_DIV (SX4, CELL_SIZE); 		// �`��(X)������3���ڊi�q�V�t�g����
static const int intSx3Per = INT_DIV (SX3, CELL_SIZE); 		// �`��(X)������3���ڊi�q�V�t�g����
static const int intSx2Per = INT_DIV (SX2, CELL_SIZE); 		// �`��(X)������3���ڊi�q�V�t�g����
static const int intSx1Per = INT_DIV (SX1, CELL_SIZE); 		// �`��(X)������1���ڊi�q�V�t�g����
static const int intSyPer = INT_DIV (SY, CELL_SIZE); 		// ��(Y)�����̓��g�H�S�̊i�q�V�t�g����


//static const int intNormPcwLen = INT_DIV (NORM_PCW_LEN, CELL_SIZE); 				// PCW����
//static const int intChirpLsLen = INT_DIV (CHIRP_3RD_LS_LEN, CELL_SIZE); 				// �`���[�vLSPCW����
//static const int intLspcwLen = INT_DIV (LSPCW_LEN, CELL_SIZE); 					// LSPCW����
//static const int intPcwLen = 2 * (intNormPcwLen + intChirpLsLen) + intLspcwLen; 	// �SPCW����
static const int intPcwWid = PCW_WID; 											// ��������PCW����
static const int intPcwStartX = intRadius; 										// �~���z�u�̊J�nX����
static const int intPcwStartY = intRadius + INT_DIV (PCW_SiSLAB_TERMINATION_LEN, CELL_SIZE);	// �~���z�u�̊J�nY����


/*****************************************************************************/
// �t�H�g�j�b�N�������g�H�ڑ��p �א����g�H
/*****************************************************************************/
static const double dblWireLen1 = WIRE_LEN1 * 1e-9; 		// ����
static const double dblWireLen2 = WIRE_LEN2 * 1e-9; 		// �o��

static const int intWireLen1 = INT_DIV(WIRE_LEN1, CELL_SIZE); 		// ����
static const int intWireLen2 = INT_DIV(WIRE_LEN2, CELL_SIZE); 		// �o��
int intWirePer2, intWirePer3; 		// �o��
static const int intWireWid_2 = intPitchY + INT_DIV(WIRE_WID_OFFSET, CELL_SIZE) + INT_DIV(SY, CELL_SIZE); 		// �א����̔����̃Z������
static const int intWireWid_2_Out = intPitchY + INT_DIV(WIRE_WID_OFFSET_OUT, CELL_SIZE) + INT_DIV(SY, CELL_SIZE);
static const int intWireWid = intWireWid_2 * 2; 		// �א�����

//#define taper_x		INT_DIV (0.0e-6, dx)		//�e�[�p����(���g�H����) Ex: (int)((0.44e-6*1e10)/(dx*1e10))
//#define taper_y		INT_DIV (0.0e-6, dx)		//�e�[�p����(������) Ex: (int)((0.24e-6*1e10)/(dy*1e10))
//#define taper		INT_DIV (0.0e-6, dx)		//�e�[�p����(���g�H�E�������ɓ������̏ꍇ) Ex: (int)((0.44e-6*1e10)/(dx*1e10))


/*****************************************************************************/
// �ޗ��̋��ܗ��ƗU�d��
/*****************************************************************************/
static const double nw = 4.5; //�ʎq���ˑw�̋��ܗ�
static const double nsch1 = 3.265; //SCH1�w�̋��ܗ� 3.265:�o���h�M���b�v�g��1100nm
static const double nsch2 = 3.292; //SCH2�w�̋��ܗ� 3.292:�o���h�M���b�v�g��1150nm
static const double nsch3 = 3.32; //SCH3�w�̋��ܗ� 3.32:�o���h�M���b�v�g��1200nm
static const double np = 3.17; //�x���̋��ܗ�
//static const double n_core = 3.5; 		// �R�A�̋��ܗ�(CORE)
static const double n_core = 3.43; 		// �R�A�̋��ܗ�(CORE)
#if PCW_Air_Or_SiO2
static const double n_clad = 1.0; 		// �N���b�h���ܗ�(�^��)
#else
static const double n_clad = 1.444; 		// �N���b�h���ܗ�(SiO2)
#endif
static const double epsilon0 = EPSILON0; 					// �^���̗U�d��
static const double epsilon1 = EPSILON0 * SQ(n_core); 		// �R�A�̗U�d��
static const double epsilon2 = EPSILON0 * SQ(n_clad); 		// �N���b�h�̗U�d��


/*****************************************************************************/
// ���U�֐�//1000
/*****************************************************************************/
#if _FDTD
#if _EXITATION_FUNC
#if _PROGRAM_TEST
#if PCW_Air_Or_SiO2
//char *dir_name[] = {"1525","1545"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1525","1545"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1525","1545","1565","1585"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1585","1588","1582","1555"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1582","1565","1575"}; 		// ���U�֐��̔g�� [nm]
char *dir_name[] = {"1565","1575"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1565"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1582"}; 		// ���U�֐��̔g�� [nm]
#else
#if PCW_S1S3_Shift
//char *dir_name[] = {"1485","1525"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1575","1565","1555","1545","1535","1525"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1565","1568","1571","1574"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1559","1556","1553","1550"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1547","1544","1541","1538"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1580","1570","1560","1550","1540"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1573","1575","1567"}; 		// ���U�֐��̔g�� [nm]
char *dir_name[] = {"1570"}; 		// ���U�֐��̔g�� [nm]
#else
//char *dir_name[] = {"1530", "1535"}; 		// ���U�֐��̔g�� [nm]
//char *dir_name[] = {"1545", "1550"}; 		// ���U�֐��̔g�� [nm]
char *dir_name[] = {"1580"}; 		// ���U�֐��̔g�� [nm]���v�Z��
#endif
#endif
#else
char *dir_name[] = {"1550"}; 		// ���U�֐��̔g�� [nm]
#endif
#else
char *dir_name[] = {"1550"}; 		// ���U�֐��̔g�� [nm]
#endif
#else
char *dir_name[] = {"1580"}; 		// ���U�֐��̔g�� [nm]���f���o����
#endif

static const double delta_omega = 0.05; 						// ���S���g���ŋK�i�������l�S��
static const int Npeak = 500; 								// �s�[�N�X�e�b�v��

// ���U�_�̍��W����
static const int ex_y_st = YMAX_ALL -24 ; 	//��-intWireWid_2	// ���g�H�f�ʎn�Z����(��) �����͋��Ԃ̒��ԃZ�����W���瓱�g�H��(����1/2�l�ɂȂ��Ă���)�������Ă��遙��
static const int ex_y_ed = YMAX_ALL; 					// ���g�H�f�ʏI�Z����(��)
static const int ex_z_st = ZMAX_ALL - intSlabHeigPer; 	// ���g�H�f�ʎn�Z����(�c) �����͋��Ԃ̒��ԃZ�����W���瓱�g�H��(��1/2�l)�������Ă���
static const int ex_z_ed = ZMAX_ALL; 			// ���g�H�f�ʏI�Z����(�c)
//static const int ex_y_st = YMAX_ALL - 26; 		// ���g�H�f�ʎn�Z����(��) �����͋��Ԃ̒��ԃZ�����W���瓱�g�H��(����1/2�l�ɂȂ��Ă���)�������Ă���
//static const int ex_y_ed = YMAX_ALL; 			// ���g�H�f�ʏI�Z����(��)
//static const int ex_z_st = ZMAX_ALL - 7; 		// ���g�H�f�ʎn�Z����(�c) �����͋��Ԃ̒��ԃZ�����W���瓱�g�H��(��1/2�l)�������Ă���
//static const int ex_z_ed = ZMAX_ALL; 			// ���g�H�f�ʏI�Z����(�c)



static const int y_zx = YMAX_ALL - 1;
/*****************************************************************************/
// �ϑ��_�Ɨ��U�_�̐ݒ� (�v�Z�덷���h�����߂Ɍ��グ���Čv�Z)
/*****************************************************************************/

// �����O�̗��U�_�C�ϑ��� (�������̓v���O�������ŏ���)
static const int intExctLen = INT_DIV (EXCT_LEN, CELL_SIZE); //���U�_ �����O�̗��U�_��
static const int intObseLen1 = INT_DIV (OBSE_LEN1, CELL_SIZE); //�ϑ��_1
//static const int intObseLen5 = XMAX_ALL - INT_DIV (OBSE_LEN5, CELL_SIZE); //���ϑ��_3 //�ϑ��_5�̊ԈႦ�ł́H�H
//�����͑����@�\���ĂȂ�

// �����O�̗��U�_�C�ϑ��� (�������̓v���O�������ŏ���)

/****************************** �ϑ��ʂ̏C��(2013/8/20) ******************************/
//static const int intObseWid = 2 * intWireWid_2; // �ϑ��ʂ̕�(Y����)
static const int intObseWid = 2 * intPitchY; // �ϑ��ʂ̕�(Y����)
static const int intObseHeig = 2 * intSlabHeigPer; // �ϑ��ʂ̍���(Z����)
/****************************** �ϑ��ʂ̏C��(2013/8/20) ******************************/

// �|�C���e�B���O�p���[�̍ő��l�C�ŏ��l�̌v�Z����
static const int intObseInter = INT_DIV (OBSE_INTER, CELL_SIZE); 	// �ϑ�����

// �������̗��U�_�C�ϑ���
int intExctPortNum;	// ���U�_�������v�Z�@�̔ԍ�
int intObseInPortNum;	// ���� �ϑ��_�������v�Z�@�̔ԍ�
int intObseOutPortNum;	// �o�� �ϑ��_�������v�Z�@�̔ԍ�
int intObseCenPortNum;	// ���� �ϑ��_�������v�Z�@�̔ԍ�

int intExctLenPart;	// ���U�_
int intObseLenPart1; // ���� �ϑ��_(�n�_)
int intObseLenPart2; // ���� �ϑ��_(�I�_)
int intObseLenPart3; // ���� �ϑ��_(���_)
int intObseLenPart4; // �o�� �ϑ��_(�n�_)
int intObseLenPart5; // �o�� �ϑ��_(�I�_)
int intObseLenPart6; // �o�� �ϑ��_(���_)
int intObseLenPart7; // ���S �ϑ��_(���_)
//int intObseLenPart8; // ���S �ϑ��_(���_)
//int intObseLenPart9; // ���S �ϑ��_(���_)

// ���E�̊ϑ���
int intObseLenPartHz1; 		// ���� �ϑ��_
int intObseLenPartHz5; 		// �o�� �ϑ��_

// �|�C���e�B���O�p���[�̍ő��l�̎Z�o�p
//static const int intObseLenPart2 = (int)((3.84e-6*1e10)/(dx*1e10)); //intObseLenPart1�����i�q�萔1��������
//static const int intObseLenPart2 = (int) INT_DIV ((2.7e-6 + dblPitch*2), dx); //intObseLenPart1�����i�q�萔2��������
//static const int intObseLenPart5 = (int)((1.32e-6*1e10)/(dx*1e10)); //intObseLenPart4�����i�q�萔1��������
//static const int intObseLenPart5 = (int)((1.3e-6*1e10 + dblPitch*2*1e10)/(dx*1e10)); //intObseLenPart4�����i�q�萔2�������� NODE 6
//static const int intObseLenPart5 = (int) INT_DIV ((9.3e-6 + dblPitch*2), dx); //intObseLenPart4�����i�q�萔2�������� NODE 2

static const int WGlength = (int)((2.00e-6*1e10)/(dx*1e10)); //���`���g�H��(���ˑ��D�Z�����ɕϊ�)
