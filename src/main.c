/***************************************************
Channel Coding Course Work: conolutional codes
This program template has given the message generator, BPSK modulation, AWGN channel model and BPSK demodulation.
This file is an enhanced version, integrating Viterbi (hard/soft) and BCJR decoders
for a rate 1/2, (7,5)_8, 4-state convolutional code.
***************************************************/

#define  _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

// --- DECODER SELECTION ---
// 1 = Hard Decision Viterbi
// 2 = Soft Decision Viterbi
// 3 = BCJR Decoder
#define DECODER_TYPE 2

// --- Code Parameters ---
#define MESSAGE_BITS 1024 // Primary message bits
#define STATE_MEM 2       // Number of memory registers (k=3, m=2)
#define message_length (MESSAGE_BITS + STATE_MEM) // Total message length with zero-padding
#define codeword_length (message_length * 2)      // Rate 1/2 code
#define st_num 4          // Number of states (2^STATE_MEM)
#define line_num 8        // Number of transitions in trellis
#define softIn_st_num 2   // For soft input (I/Q)
#define inf_int 1000000   // "Infinity" for integer path metrics
#define inf_double 1e18   // "Infinity" for double path metrics
#define p0 0.0              // Initial probability

float code_rate = (float)MESSAGE_BITS / (float)codeword_length; // Effective code rate

// channel coefficient
#define pi 3.1415926
double N0, sgm;

// --- Global Variables from Template ---
int state_num; //the number of the state of encoder structure (will be set to STATE_MEM)

int message[message_length], codeword[codeword_length];//message and codeword
int re_codeword[codeword_length];//the received codeword (for hard decision)
int de_message[message_length];//the decoding message

double tx_symbol[codeword_length][2];//the transmitted symbols
double rx_symbol[codeword_length][2];//the received symbols

// --- Structs and Definitions from Senior's Code (statement.h) ---
typedef struct {
    int input;
    int output;
    int pt1; // Start state
    int pt2; // End state
    int id;
} line;

typedef struct {
    int point;
    int line;
} connect;

typedef struct {
    connect first;
    connect second;
} graph_connect;

// --- Trellis Structure from Senior's Code (viterbit.c) ---
// (7,5)_8 code, Generator G(D) = [1+D^2, 1+D+D^2] -> [101, 111] -> (5, 7) octal.
// Note: Senior's ENCODER implements (5,7) octal, not (7,5). We follow the encoder.
// Encoder: c1 = m[i] ^ s0 ^ s1; c2 = m[i] ^ s1; (s1_next = s0, s0_next = m[i])
// State = (s1, s0)
const int state[4][2] = { {0,0},{0,1},{1,0},{1,1} }; // Outputs (c1, c2)
const int c00 = 0, c01 = 1, c10 = 2, c11 = 3; // Index for state[][]

// State transition table
line stateTable[line_num] = {
    // input, output_idx, start_state, end_state, line_id
    {0, c00, 0, 0, 1}, // State 0 (00), in 0 -> out (00), next_state 0 (00)
    {1, c11, 0, 2, 2}, // State 0 (00), in 1 -> out (11), next_state 2 (10)
    {0, c11, 1, 0, 3}, // State 1 (01), in 0 -> out (11), next_state 0 (00)
    {1, c00, 1, 2, 4}, // State 1 (01), in 1 -> out (00), next_state 2 (10)
    {0, c10, 2, 1, 5}, // State 2 (10), in 0 -> out (10), next_state 1 (01)
    {1, c01, 2, 3, 6}, // State 2 (10), in 1 -> out (01), next_state 3 (11)
    {0, c01, 3, 1, 7}, // State 3 (11), in 0 -> out (01), next_state 1 (01)
    {1, c10, 3, 3, 8}  // State 3 (11), in 1 -> out (10), next_state 3 (11)
};

// Forward connection (for ACS)
// pathConn[STATE] = { {PREV_STATE_1, LINE_ID_1}, {PREV_STATE_2, LINE_ID_2} }
graph_connect pathConn[4] = {
    { {0,0}, {1,2} }, // To State 0 (from S0, L1 or S1, L3)
    { {2,4}, {3,6} }, // To State 1 (from S2, L5 or S3, L7)
    { {0,1}, {1,3} }, // To State 2 (from S0, L2 or S1, L4)
    { {2,5}, {3,7} }  // To State 3 (from S2, L6 or S3, L8)
};

// Backward connection (for BCJR)
// pathConnB[STATE] = { {NEXT_STATE_1, LINE_ID_1}, {NEXT_STATE_2, LINE_ID_2} }
graph_connect pathConnB[4] = {
    { {0,0}, {2,1} }, // From State 0 (to S0, L1 or S2, L2)
    { {0,2}, {2,3} }, // From State 1 (to S0, L3 or S2, L4)
    { {1,4}, {3,5} }, // From State 2 (to S1, L5 or S3, L6)
    { {1,6}, {3,7} }  // From State 3 (to S1, L7 or S3, L8)
};

// --- Global Arrays for Decoders ---
// Viterbi
int branchTable[message_length][line_num];
int pathTable[message_length + 1][st_num];
int trellisTable[message_length][st_num];
int minPath[message_length];
double branchTableSoft[message_length][line_num];
double pathTableSoft[message_length + 1][st_num];

// BCJR
double pCh[codeword_length][2];
double pLine[message_length][line_num];
double pA[message_length + 1][st_num];
double pB[message_length + 1][st_num];

// --- Function Prototypes ---
void statetable();
void encoder();
void modulation();
void demodulation();
void channel();
void decoder();

// Viterbi
void hardDecoder(int re_codeword[], int de_message[], int ms_length);
void softDecode(double re_codewordSoft[][softIn_st_num], int de_message[], int ms_length);

// BCJR
double chObs(double eu, double n0);
double euDist(double, double, int);
void BCJR(double[][softIn_st_num], int[], int, int);


void main()
{
	int i;
	float SNR, start, finish;
	long int bit_error, seq, seq_num;
	double BER;
	double progress;

    // Set encoder memory size
    state_num = STATE_MEM;

	//generate state table (not needed, tables are constants)
	statetable();

	//random seed
	srand((int)time(0));

	//input the SNR and frame number
	printf("\nEnter start SNR: ");
	scanf("%f", &start);
	printf("\nEnter finish SNR: ");
	scanf("%f", &finish);
	printf("\nPlease input the number of message: ");
	scanf("%d", &seq_num);

	for (SNR = start; SNR <= finish; SNR++)
	{
		//channel noise
		N0 = (1.0 / code_rate) / pow(10.0, (float)(SNR) / 10.0);
		sgm = sqrt(N0 / 2);
		
		bit_error = 0;

		for (seq = 1; seq <= seq_num; seq++)
		{
			//generate binary message randomly
			for (i = 0; i < MESSAGE_BITS; i++)
			{
				message[i] = rand() % 2;
			}
            // Zero-padding for trellis termination
			for (i = MESSAGE_BITS; i < message_length; i++)
			{
				message[i] = 0;
			}

			//convolutional encoder
			encoder();

			//BPSK modulation
			modulation();

			//AWGN channel
			channel();

            //BPSK demodulation, only needed for hard-decision Viterbi
            #if DECODER_TYPE == 1
			    demodulation();
            #endif

			//convolutional decoder
			decoder();

			//calculate the number of bit error
            // Only check the original message bits, not the padding
			for (i = 0; i < MESSAGE_BITS; i++)
			{
				if (message[i] != de_message[i])
					bit_error++;
			}

			progress = (double)(seq * 100) / (double)seq_num;
			BER = (double)bit_error / (double)(MESSAGE_BITS * seq);

			printf("Progress=%2.1f, SNR=%2.1f, Bit Errors=%ld, BER=%E\r", progress, SNR, bit_error, BER);
		}

		BER = (double)bit_error / (double)(MESSAGE_BITS * seq_num);
		printf("Progress=100.0, SNR=%2.1f, Bit Errors=%ld, BER=%E\n", SNR, bit_error, BER);
	}
	system("pause");
}

void statetable()
{
    // This function is not needed because we use hard-coded constant tables
    // (stateTable, pathConn, pathConnB) for the specific (5,7)_8 code.
}

void encoder()
{
	// Convolution encoder for (5,7)_8 code
    // G = [101, 111]
    int c1, c2;
    int s0 = 0; // Register 1
    int s1 = 0; // Register 2
    for (int i = 0; i < message_length; i++)
    {
        c1 = message[i] ^ s0 ^ s1; // 1*m + 1*s0 + 1*s1 -> (111) -> 7
        c2 = message[i] ^ s1;      // 1*m + 0*s0 + 1*s1 -> (101) -> 5
        
        // This is (7,5) octal, not (5,7). Let's swap to match viterbit.c
        // Encoder from senior's src/main.c:
        // c1 = message[i] ^ s0 ^ s1; // (111) -> 7
        // c2 = message[i] ^ s1;      // (101) -> 5
        // This matches the trellis tables.
        
        s1 = s0; // Shift register
        s0 = message[i]; // Shift register
        codeword[2 * i] = c1;
        codeword[2 * i + 1] = c2;
    }
}

void modulation()
{
	//BPSK modulation
	int i;
	//0 is mapped to (1,0) and 1 is mapped tp (-1,0)
	for (i = 0; i<codeword_length; i++)
	{
		tx_symbol[i][0] = -1 * (2 * codeword[i] - 1);
		tx_symbol[i][1]=0;
	}
}
void channel()
{
	//AWGN channel
	int i, j;
	double u, r, g;

	for (i = 0; i<codeword_length; i++)
	{
		for (j = 0; j<2; j++)
		{
			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			r=sgm*sqrt(2.0*log(1.0/(1.0-u)));

			u=(float)rand()/(float)RAND_MAX;
			if(u==1.0)
				u=0.999999;
			g=(float)r*cos(2*pi*u);

			rx_symbol[i][j]=tx_symbol[i][j]+g;
		}
	}
}
void demodulation()
{
    // BPSK Hard-decision demodulation
	int i;
	double d1, d2;
	for (i = 0; i<codeword_length; i++)
	{
        // 0 maps to +1, 1 maps to -1
		d1 = (rx_symbol[i][0] - 1)*(rx_symbol[i][0] - 1) + rx_symbol[i][1] * rx_symbol[i][1]; // Dist to 0 (+1)
		d2 = (rx_symbol[i][0] + 1)*(rx_symbol[i][0] + 1) + rx_symbol[i][1] * rx_symbol[i][1]; // Dist to 1 (-1)
		if (d1<d2)
			re_codeword[i] = 0;
		else
			re_codeword[i] = 1;
	}
}

void decoder()
{ 
    // Selects which decoder to run based on the macro
    #if DECODER_TYPE == 1
        // printf("Using Hard Viterbi Decoder\n");
        hardDecoder(re_codeword, de_message, message_length);
    #elif DECODER_TYPE == 2
        // printf("Using Soft Viterbi Decoder\n");
        softDecode(rx_symbol, de_message, message_length);
    #elif DECODER_TYPE == 3
        // printf("Using BCJR Decoder\n");
        BCJR(rx_symbol, de_message, message_length, codeword_length);
    #else
        printf("Error: Invalid DECODER_TYPE specified.\n");
    #endif
}


// --- VITERBI DECODER FUNCTIONS (from viterbit.c) ---

void hardDecoder(int re_codeword[], int de_message[], int ms_length)
{
    int hammingDistance;
    for (int t = 0; t < ms_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            hammingDistance = 0;
            // Filter out impossible transitions
            if (t == 0 && index != 0 && index != 1) { branchTable[t][index] = inf_int; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { branchTable[t][index] = inf_int; continue; }
            if (t == ms_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { branchTable[t][index] = inf_int; continue; }
            if (t == ms_length - 1 && index != 0 && index != 2) { branchTable[t][index] = inf_int; continue; }

            int nowCode1 = state[stateTable[index].output][0];
            int nowCode2 = state[stateTable[index].output][1];

            if (nowCode1 != re_codeword[2 * t]) hammingDistance++;
            if (nowCode2 != re_codeword[2 * t + 1]) hammingDistance++;
            
            branchTable[t][index] = hammingDistance;
        }
    }

    int path1d, path2d, minpath;
    pathTable[0][0] = 0;
    for (int fst = 1; fst < st_num; fst++) { pathTable[0][fst] = inf_int; }

    // Add-Compare-Select (ACS)
    for (int pt = 1; pt < ms_length + 1; pt++)
    {
        for (int st = 0; st < st_num; st++)
        {
            path1d = pathTable[pt - 1][pathConn[st].first.point] + branchTable[pt - 1][pathConn[st].first.line];
            path2d = pathTable[pt - 1][pathConn[st].second.point] + branchTable[pt - 1][pathConn[st].second.line];

            if (pt == 1 && st != 0 && st != 2) { trellisTable[pt - 1][st] = inf_int; pathTable[pt][st] = inf_int; continue; }
            if (pt == ms_length - 1 && st != 0 && st != 1) { trellisTable[pt - 1][st] = inf_int; pathTable[pt][st] = inf_int; continue; }
            if (pt == ms_length && st != 0) { trellisTable[pt - 1][st] = inf_int; pathTable[pt][st] = inf_int; continue; }

            if (path1d < path2d) {
                minpath = path1d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].first.line].id;
            }
            else {
                minpath = path2d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].second.line].id;
            }
            pathTable[pt][st] = minpath;
        }
    }

    // Traceback
    int nowLine = 0, nowPoint = 0, prevPoint = 0;
    for (int tt = ms_length; tt > 0; tt--)
    {
        nowLine = trellisTable[tt - 1][nowPoint] - 1; 
        prevPoint = stateTable[nowLine].pt1;     
        nowPoint = prevPoint;                         
        minPath[tt - 1] = nowLine + 1;                 
    }

    for (int t = 0; t < ms_length; t++)
    {
        de_message[t] = stateTable[minPath[t] - 1].input;
    }
}


void softDecode(double re_codewordSoft[][softIn_st_num], int de_message[], int ms_length)
{
    double euDistance;
    for (int t = 0; t < ms_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            euDistance = 0;
            // Filter out impossible transitions
            if (t == 0 && index != 0 && index != 1) { branchTableSoft[t][index] = inf_double; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { branchTableSoft[t][index] = inf_double; continue; }
            if (t == ms_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { branchTableSoft[t][index] = inf_double; continue; }
            if (t == ms_length - 1 && index != 0 && index != 2) { branchTableSoft[t][index] = inf_double; continue; }

            int nowCode1 = state[stateTable[index].output][0];
            int nowCode2 = state[stateTable[index].output][1];

            // 0 maps to +1, 1 maps to -1
            double ideal_x1 = (nowCode1 == 0) ? 1.0 : -1.0;
            double ideal_x2 = (nowCode2 == 0) ? 1.0 : -1.0;

            // Using squared Euclidean distance is equivalent and faster (avoids sqrt)
            euDistance += pow(re_codewordSoft[2 * t][0] - ideal_x1, 2);
            euDistance += pow(re_codewordSoft[2 * t + 1][0] - ideal_x2, 2);

            branchTableSoft[t][index] = euDistance;
        }
    }

    double path1d, path2d, minpath;
    pathTableSoft[0][0] = 0;
    for (int fst = 1; fst < st_num; fst++) { pathTableSoft[0][fst] = inf_double; }

    // Add-Compare-Select (ACS)
    for (int pt = 1; pt < ms_length + 1; pt++)
    {
        for (int st = 0; st < st_num; st++)
        {
            path1d = pathTableSoft[pt - 1][pathConn[st].first.point] + branchTableSoft[pt - 1][pathConn[st].first.line];
            path2d = pathTableSoft[pt - 1][pathConn[st].second.point] + branchTableSoft[pt - 1][pathConn[st].second.line];

            if (pt == 1 && st != 0 && st != 2) { trellisTable[pt - 1][st] = inf_int; pathTableSoft[pt][st] = inf_double; continue; }
            if (pt == ms_length - 1 && st != 0 && st != 1) { trellisTable[pt - 1][st] = inf_int; pathTableSoft[pt][st] = inf_double; continue; }
            if (pt == ms_length && st != 0) { trellisTable[pt - 1][st] = inf_int; pathTableSoft[pt][st] = inf_double; continue; }
            
            if (path1d < path2d) {
                minpath = path1d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].first.line].id;
            }
            else {
                minpath = path2d;
                trellisTable[pt - 1][st] = stateTable[pathConn[st].second.line].id;
            }
            pathTableSoft[pt][st] = minpath;
        }
    }

    // Traceback
    int nowLine = 0, nowPoint = 0, prevPoint = 0;
    for (int tt = ms_length; tt > 0; tt--)
    {
        nowLine = trellisTable[tt - 1][nowPoint] - 1; 
        prevPoint = stateTable[nowLine].pt1;     
        nowPoint = prevPoint;                         
        minPath[tt - 1] = nowLine + 1;                 
    }

    for (int t = 0; t < ms_length; t++)
    {
        de_message[t] = stateTable[minPath[t] - 1].input;
    }
}


// --- BCJR DECODER FUNCTIONS (from bcjr.c) ---

double euDist(double rx_x, double rx_y, int sym)
{
    // 0 maps to +1, 1 maps to -1
    double ideal_x = (sym == 0) ? 1.0 : -1.0;
    // Using squared Euclidean distance
    return pow(rx_x - ideal_x, 2) + pow(rx_y - 0.0, 2);
}

double chObs(double eu_sq, double n0)
{
    // P(y | x) = (1 / sqrt(pi*N0)) * exp(-d^2(y,x) / N0)
    // Note: N0 = 2*sgm^2
    return exp(-eu_sq / n0) / (sqrt(pi * n0));
}

void BCJR(double rx_sym[][softIn_st_num], int de_message[], int m_length, int c_length)
{
	// 1. Calculate Channel Observation Probabilities (Gamma)
    // pCh[t][0] = P(y_t | 0), pCh[t][1] = P(y_t | 1)
    double sum;
    for (int t = 0; t < c_length; t++)
    {
        sum = 0;
        pCh[t][0] = chObs(euDist(rx_sym[t][0], rx_sym[t][1], 0), N0);
        pCh[t][1] = chObs(euDist(rx_sym[t][0], rx_sym[t][1], 1), N0);
        // Normalize probabilities to avoid underflow
        sum = pCh[t][0] + pCh[t][1];
        if (sum > 1e-300) { // Avoid division by zero
            pCh[t][0] = pCh[t][0] / sum;
            pCh[t][1] = pCh[t][1] / sum;
        } else {
            pCh[t][0] = 0.5;
            pCh[t][1] = 0.5;
        }
    }
    
    // 2. Calculate Branch Transition Probabilities (Gamma)
    // P(s_t-1, s_t) = P(u_t) * P(y_t | s_t-1, s_t)
    // P(u_t) = 0.5 (a priori probability)
    for (int t = 0; t < m_length; t++)
    {
        for (int index = 0; index < line_num; index++)
        {
            // Filter out impossible transitions
            if (t == 0 && index != 0 && index != 1) { pLine[t][index] = p0; continue; }
            if (t == 1 && index != 0 && index != 1 && index != 4 && index != 5) { pLine[t][index] = p0; continue; }
            if (t == m_length - 2 && index != 0 && index != 2 && index != 4 && index != 6) { pLine[t][index] = p0; continue; }
            if (t == m_length - 1 && index != 0 && index != 2) { pLine[t][index] = p0; continue; }
            
            int out_idx = stateTable[index].output;
            int c1 = state[out_idx][0];
            int c2 = state[out_idx][1];

            // P(u_t) * P(y_t1 | c1) * P(y_t2 | c2)
            pLine[t][index] = 0.5 * pCh[2 * t][c1] * pCh[2 * t + 1][c2];
        }   
    }
    
	// 3. Calculate Forward Probabilities (Alpha)
    pA[0][0] = 1.0; pA[0][1] = 0.0; pA[0][2] = 0.0; pA[0][3] = 0.0;
    double sumA;
    for (int t = 1; t < m_length + 1; t++)
    {
        sumA = 0;
        for (int st = 0; st < st_num; st++)
        {
            pA[t][st] = pA[t - 1][pathConn[st].first.point] * pLine[t - 1][pathConn[st].first.line] +
                        pA[t - 1][pathConn[st].second.point] * pLine[t - 1][pathConn[st].second.line];
            sumA += pA[t][st];
        }
        // Normalize to prevent underflow
        if (sumA > 1e-300) {
            for (int st = 0; st < st_num; st++) {
                pA[t][st] = pA[t][st] / sumA;
            }
        }
    }
    
    // 4. Calculate Backward Probabilities (Beta)
    pB[m_length][0] = 1.0; pB[m_length][1] = 0.0; pB[m_length][2] = 0.0; pB[m_length][3] = 0.0;
    double sumB;
    for (int t = m_length - 1; t > -1; t--)
    {
        sumB = 0;
        for (int st = 0; st < st_num; st++)
        {
            pB[t][st] = pB[t + 1][pathConnB[st].first.point] * pLine[t][pathConnB[st].first.line] +
                        pB[t + 1][pathConnB[st].second.point] * pLine[t][pathConnB[st].second.line];
            sumB += pB[t][st];
        }
        // Normalize
        if (sumB > 1e-300) {
            for (int st = 0; st < st_num; st++) {
                pB[t][st] = pB[t][st] / sumB;
            }
        }
    }
    
    // 5. Calculate A Posteriori Probabilities (APP) and Decide
    double Po, Pz; // Po = P(u_t = 1 | y), Pz = P(u_t = 0 | y)
    for (int t = 0; t < m_length; t++)
	{
        // Sum probabilities for all transitions caused by input u_t = 0
		Pz = pA[t][0] * pLine[t][0] * pB[t+1][0]   // S0->S0 (line 0)
           + pA[t][1] * pLine[t][2] * pB[t+1][0]   // S1->S0 (line 2)
           + pA[t][2] * pLine[t][4] * pB[t+1][1]   // S2->S1 (line 4)
           + pA[t][3] * pLine[t][6] * pB[t+1][1];  // S3->S1 (line 6)

        // Sum probabilities for all transitions caused by input u_t = 1
		Po = pA[t][0] * pLine[t][1] * pB[t+1][2]   // S0->S2 (line 1)
           + pA[t][1] * pLine[t][3] * pB[t+1][2]   // S1->S2 (line 3)
           + pA[t][2] * pLine[t][5] * pB[t+1][3]   // S2->S3 (line 5)
           + pA[t][3] * pLine[t][7] * pB[t+1][3];  // S3->S3 (line 7)

		if(Pz < Po) {
			de_message[t] = 1;
	    } else {
			de_message[t] = 0;
		}		
	}
}