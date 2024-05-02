/**
 *******************************************************************************
 * @version 0.1
 * @date 2021-08-22
 * @copyright Copyright © 2021 by University of Luxembourg.
 * @author Developed at SnT APSIA by: Hao Cheng.
 *******************************************************************************
 */

#include "action.h"
#include "utils.h"
#include <string.h>

// the function to measure CPU cycles 
extern uint64_t read_tsc();

// MACROs for benchmarking 

#define ITER_L 100000
#define ITER_M 10000
#define ITER_S 1000

#define LOAD_CACHE(X, ITER) for (i = 0; i < (ITER); i++) (X)

#define MEASURE_TIME(X, ITER)                         \
  start_cycles = read_tsc();                          \
  for (i = 0; i < (ITER); i++) (X);                   \
  end_cycles = read_tsc();                            \
  diff_cycles = (end_cycles-start_cycles)/(ITER)

static void sk_print(char *c, uint8_t *sk)
{
  int i;

  printf("%s", c);
  printf("%d", (int)((2*(sk[0]&1)-1)*(sk[0]>>1)));

  for(i = 1; i < N; i++)
  {
    printf(":%d", (int)((2*(sk[i]&0x1)-1)*(sk[i]>>1)));
  };
  printf("\n");
};

void test_action()
{
  // test vectors for OAYT-style implementation
  // have been verified by [CCC+19] code
  // - SKA = -2:-3:1:3:-2:0:-2:3:3:-1:3:0:-1:1:-4:1:4:1:2:0:0:-3:4:4:-2:0:0:-2:5:-1:5:-1:4:4:5:3:6:5:1:-5:-7:-3:6:-7:-1:-7:-2:1:-3:-2:2:1:-3:7:-10:2:-3:8:5:2:3:7:4:-5:-1:-3:-1:2:2:0:-1:0:1:1
  // - SKB = 1:-1:-3:-1:0:3:3:-1:-2:3:2:0:-4:-3:-3:4:-3:-3:2:-3:3:0:2:-4:-2:-3:3:4:0:3:5:-2:2:-6:-6:-2:1:4:1:1:4:1:-5:-1:-7:-5:-5:-5:7:6:8:4:5:9:10:1:-5:10:-2:8:-4:0:-1:5:3:7:3:5:2:0:1:0:0:1
  // - PKA : 08A0E2B7EC266A371ECC73F2661217AA507A2986EDC109E9F01C0E47D62BAC0069F5132AAC7C1CAB2F91DAF8A827729D70AFB5EA93373AB0D7F5B69014B6F7ED
  // - PKB : 3184BCFFDB31F8B626CD0CFD99D30D3E1CD83F7F46BB50870A433F3F39478ACC830D64647930E1A093FBA8CD221576AD1C33DCF3279F8C84A7770D9E569423D9
  // - SS  : 5647C828E4C8112F1698C130B689946408674D3E264797F5D3C18D95DC579BB8E69EAB7018C7B41EA58B38F055EF72B46636EF9850E7CA6053670C8A16C918C3

  // uncomment the following two lines if you want to use test vectors for OAYT-style implementation 
  // uint8_t sk_a[N] = {4, 6, 3, 7, 4, 1, 4, 7, 7, 2, 7, 1, 2, 3, 8, 3, 9, 3, 5, 1, 1, 6, 9, 9, 4, 1, 1, 4, 11, 2, 11, 2, 9, 9, 11, 7, 13, 11, 3, 10, 14, 6, 13, 14, 2, 14, 4, 3, 6, 4, 5, 3, 6, 15, 20, 5, 6, 17, 11, 5, 7, 15, 9, 10, 2, 6, 2, 5, 5, 1, 2, 1, 3, 3, };
  // uint8_t sk_b[N] = {3, 2, 6, 2, 1, 7, 7, 2, 4, 7, 5, 1, 8, 6, 6, 9, 6, 6, 5, 6, 7, 1, 5, 8, 4, 6, 7, 9, 1, 7, 11, 4, 5, 12, 12, 4, 3, 9, 3, 3, 9, 3, 10, 2, 14, 10, 10, 10, 15, 13, 17, 9, 11, 19, 21, 3, 10, 21, 4, 17, 8, 1, 2, 11, 7, 15, 7, 11, 5, 1, 3, 1, 1, 3, };
  
  // test vectors for Dummy-Free-style implementation
  // have been verified by [CCC+19] code
  // - SKA = 7:-5:3:1:3:5:-5:1:-3:1:1:3:-3:5:1:-5:-8:0:-2:8:-8:4:8:11:-1:5:-5:-3:7:-3:-3:3:-7:-5:5:-7:-5:3:-5:-9:-1:-11:11:13:7:-11:13:-1:13:-9:-5:13:-1:11:11:-5:-7:9:7:5:1:1:-5:-11:5:-7:7:-7:-5:-1:7:7:-7:7
  // - SKB = -5:3:-7:-5:-5:5:3:5:-1:-5:-1:5:7:7:3:7:-8:6:2:-8:-8:-6:8:-11:-9:-1:-9:-1:11:-5:3:-5:1:-1:5:3:3:-7:-9:-9:1:5:5:-5:7:-7:-13:13:-9:-5:3:11:-7:3:-7:-1:11:-5:3:7:-9:9:-11:-3:-7:11:5:-9:-1:5:5:-3:7:3
  // - PKA = 184F578A202B1FBFC0257243DBE9810569D7858E5D8239A2D682258F6579F00DF33F98CB06D14A9F8752E1BCADD527D7A9B619E2E5B8AFCA38B5FE031731C10E
  // - PKB = 61A8176234DACCFCE6C618E7BBECA46BD38DD5E6C67FD67D5B9E8EF0DD93478147AE17D0CF8C0C4DABB317B1C85147E2FF6869994BD9A15ADF947BB683672EC6  
  // - SS  = 09E943B1955EBE3A8FEDE1137FC17BA3859E581F53522A364670EBBA02CAC4358CE2EEEC138C13814053CEC9C68AF312B7886460CD9195599D854092BDCB92FC

  // uncomment the following two lines if you want to use test vectors for Dummy-Free-style implementation 
  // uint8_t sk_a[N] = {15, 10, 7, 3, 7, 11, 10, 3, 6, 3, 3, 7, 6, 11, 3, 10, 16, 1, 4, 17, 16, 9, 17, 23, 2, 11, 10, 6, 15, 6, 6, 7, 14, 10, 11, 14, 10, 7, 10, 18, 2, 22, 23, 27, 15, 22, 27, 2, 27, 18, 10, 27, 2, 23, 23, 10, 14, 19, 15, 11, 3, 3, 10, 22, 11, 14, 15, 14, 10, 2, 15, 15, 14, 15};
  // uint8_t sk_b[N] = {10, 7, 14, 10, 10, 11, 7, 11, 2, 10, 2, 11, 15, 15, 7, 15, 16, 13, 5, 16, 16, 12, 17, 22, 18, 2, 18, 2, 23, 10, 7, 10, 3, 2, 11, 7, 7, 14, 18, 18, 3, 11, 11, 10, 15, 14, 26, 27, 18, 10, 7, 23, 14, 7, 14, 2, 23, 10, 7, 15, 18, 19, 22, 6, 14, 23, 11, 18, 2, 11, 11, 6, 15, 7};


  uint8_t sk_a[N], sk_b[N];
  uint8_t sk_c[N], sk_d[N];
  uint8_t sk_e[N], sk_f[N];
  uint8_t sk_g[N], sk_h[N];
  __m512i vsk[N];
  htpoint vPK, vSS, vE;
  htfe_t t0, t1;
  uint32_t pk_a[HT_NWORDS], pk_b[HT_NWORDS], pk_c[HT_NWORDS], pk_d[HT_NWORDS];
  uint32_t pk_e[HT_NWORDS], pk_f[HT_NWORDS], pk_g[HT_NWORDS], pk_h[HT_NWORDS];
  uint32_t ss_a[HT_NWORDS] = { 0 }, ss_b[HT_NWORDS] = { 0 };
  uint32_t ss_c[HT_NWORDS] = { 0 }, ss_d[HT_NWORDS] = { 0 };
  uint32_t ss_e[HT_NWORDS] = { 0 }, ss_f[HT_NWORDS] = { 0 };
  uint32_t ss_g[HT_NWORDS] = { 0 }, ss_h[HT_NWORDS] = { 0 };
  uint32_t r29[HT_NWORDS], r32[HT_NWORDS];
  int i;

  puts("\n*******************************************************************");
  puts("CORRECTNESS TEST:");
  puts("-------------------------------------------------------------------");
  puts("EIGHT instances in each execution of the program.");
  puts("We can test our software by assuming four pairs of them are sharing the secret:");
  puts("A  <---------------------------->  B");
  puts("C  <---------------------------->  D");
  puts("E  <---------------------------->  F");
  puts("G  <---------------------------->  H");

  // randomly generate eight private keys (comment the following line if you are using test vectors)
  random_sk(sk_a); random_sk(sk_b); 
  random_sk(sk_c); random_sk(sk_d);
  random_sk(sk_e); random_sk(sk_f); 
  random_sk(sk_g); random_sk(sk_h);

  // form the vector of private key
  for (i = 0; i < N; i++) 
    vsk[i] = set_vector(sk_h[i], sk_g[i], sk_f[i], sk_e[i], sk_d[i], sk_c[i], sk_b[i], sk_a[i]);

  // intialize the vector of base curve
  for (i = 0; i < HT_NWORDS; i++) {
    vE.y[i] = VSET1(E[0][i]);
    vE.z[i] = VSET1(E[1][i]);
  }

  puts("\n* Privake Keys:");
  sk_print("  - SKA = ", sk_a);
  sk_print("  - SKB = ", sk_b);
  sk_print("  - SKC = ", sk_c);
  sk_print("  - SKD = ", sk_d);
  sk_print("  - SKE = ", sk_e);
  sk_print("  - SKF = ", sk_f);
  sk_print("  - SKG = ", sk_g);
  sk_print("  - SKH = ", sk_h);

  puts("\n* Public Keys:");  
  action(&vPK, vsk, &vE);

  gfp_inv_8x1w(t0, vPK.z);
  gfp_mul_8x1w(t1, vPK.y, t0);
  gfp_mont2num_8x1w(t0, t1);

  get_channel_8x1w(pk_a, t0, 0);
  get_channel_8x1w(pk_b, t0, 1);
  get_channel_8x1w(pk_c, t0, 2);
  get_channel_8x1w(pk_d, t0, 3);
  get_channel_8x1w(pk_e, t0, 4);
  get_channel_8x1w(pk_f, t0, 5);
  get_channel_8x1w(pk_g, t0, 6);
  get_channel_8x1w(pk_h, t0, 7);

  mpi_conv_29to32(r32, pk_a, HT_NWORDS, HT_NWORDS);
  mpi_print("  - A : ", r32, 16);
  mpi_conv_29to32(r32, pk_b, HT_NWORDS, HT_NWORDS);
  mpi_print("  - B : ", r32, 16);
  mpi_conv_29to32(r32, pk_c, HT_NWORDS, HT_NWORDS);
  mpi_print("  - C : ", r32, 16);
  mpi_conv_29to32(r32, pk_d, HT_NWORDS, HT_NWORDS);
  mpi_print("  - D : ", r32, 16);
  mpi_conv_29to32(r32, pk_e, HT_NWORDS, HT_NWORDS);
  mpi_print("  - E : ", r32, 16);
  mpi_conv_29to32(r32, pk_f, HT_NWORDS, HT_NWORDS);
  mpi_print("  - F : ", r32, 16);
  mpi_conv_29to32(r32, pk_g, HT_NWORDS, HT_NWORDS);
  mpi_print("  - G : ", r32, 16);
  mpi_conv_29to32(r32, pk_h, HT_NWORDS, HT_NWORDS);
  mpi_print("  - H : ", r32, 16);

  // We simply swap the position of (pk_a and pk_b), (pk_c and pk_d),
  // (pk_e and pk_f), (pk_g and pk_h) in the AVX512 vectors. 
  // So now it is (pk_g, pk_h, pk_e, pk_f, pk_c, pk_d, pk_a, pk_b)
  //              (sk_h, sk_g, sk_f, sk_e, sk_d, sk_c, sk_b, sk_a)
  for (i = 0; i < HT_NWORDS; i++) { 
    vPK.y[i] = VPERM(vPK.y[i], 0xB1);
    vPK.z[i] = VPERM(vPK.z[i], 0xB1);
  }

  puts("\n* Shared Secrets:");  
  action(&vSS, vsk, &vPK);

  gfp_inv_8x1w(t0, vSS.z);
  gfp_mul_8x1w(t1, vSS.y, t0);
  gfp_mont2num_8x1w(t0, t1);

  get_channel_8x1w(ss_a, t0, 0);
  get_channel_8x1w(ss_b, t0, 1);
  get_channel_8x1w(ss_c, t0, 2);
  get_channel_8x1w(ss_d, t0, 3);
  get_channel_8x1w(ss_e, t0, 4);
  get_channel_8x1w(ss_f, t0, 5);
  get_channel_8x1w(ss_g, t0, 6);
  get_channel_8x1w(ss_h, t0, 7);
  mpi_conv_29to32(r32, ss_a, HT_NWORDS, HT_NWORDS);
  mpi_print("  - A : ", r32, 16);
  mpi_conv_29to32(r32, ss_b, HT_NWORDS, HT_NWORDS);
  mpi_print("  - B : ", r32, 16);
  mpi_conv_29to32(r32, ss_c, HT_NWORDS, HT_NWORDS);
  mpi_print("  - C : ", r32, 16);
  mpi_conv_29to32(r32, ss_d, HT_NWORDS, HT_NWORDS);
  mpi_print("  - D : ", r32, 16);
  mpi_conv_29to32(r32, ss_e, HT_NWORDS, HT_NWORDS);
  mpi_print("  - E : ", r32, 16);
  mpi_conv_29to32(r32, ss_f, HT_NWORDS, HT_NWORDS);
  mpi_print("  - F : ", r32, 16);
  mpi_conv_29to32(r32, ss_g, HT_NWORDS, HT_NWORDS);
  mpi_print("  - G : ", r32, 16);
  mpi_conv_29to32(r32, ss_h, HT_NWORDS, HT_NWORDS);
  mpi_print("  - H : ", r32, 16);

  puts("\n* Correctness:");
  if (memcmp(ss_a, ss_b, HT_NWORDS*sizeof(uint32_t)))
        printf("Shared secret between A and B : \x1b[31mNOT EQUAL!\x1b[0m\n");
  else
    printf("Shared secret between A and B : \x1b[32mEQUAL!\x1b[0m\n");

  if (memcmp(ss_c, ss_d, HT_NWORDS*sizeof(uint32_t)))
    printf("Shared secret between C and D : \x1b[31mNOT EQUAL!\x1b[0m\n");
  else
    printf("Shared secret between C and D : \x1b[32mEQUAL!\x1b[0m\n");

  if (memcmp(ss_e, ss_f, HT_NWORDS*sizeof(uint32_t)))
        printf("Shared secret between E and F : \x1b[31mNOT EQUAL!\x1b[0m\n");
  else
    printf("Shared secret between E and F : \x1b[32mEQUAL!\x1b[0m\n");

  if (memcmp(ss_g, ss_h, HT_NWORDS*sizeof(uint32_t)))
    printf("Shared secret between G and H : \x1b[31mNOT EQUAL!\x1b[0m\n");
  else
    printf("Shared secret between G and H : \x1b[32mEQUAL!\x1b[0m\n");

  puts("*******************************************************************");
}

void test_multi_actions(int round)
{
  int i, j, wrong = 0;
  uint8_t sk_a[N], sk_b[N], sk_c[N], sk_d[N];
  uint8_t sk_e[N], sk_f[N], sk_g[N], sk_h[N];
  __m512i vsk[N];
  htpoint vPK, vSS, vE;
  htfe_t t0, t1;
  uint32_t pk_a[HT_NWORDS], pk_b[HT_NWORDS], pk_c[HT_NWORDS], pk_d[HT_NWORDS];
  uint32_t pk_e[HT_NWORDS], pk_f[HT_NWORDS], pk_g[HT_NWORDS], pk_h[HT_NWORDS];
  uint32_t ss_a[HT_NWORDS] = { 0 }, ss_b[HT_NWORDS] = { 0 };
  uint32_t ss_c[HT_NWORDS] = { 0 }, ss_d[HT_NWORDS] = { 0 };
  uint32_t ss_e[HT_NWORDS] = { 0 }, ss_f[HT_NWORDS] = { 0 };
  uint32_t ss_g[HT_NWORDS] = { 0 }, ss_h[HT_NWORDS] = { 0 };
  uint32_t r29[HT_NWORDS], r32[HT_NWORDS];

  puts("\n* Multi Tests:");

  for (j = 0; j < round; j++) {

    // randomly generate eight private keys 
    random_sk(sk_a); random_sk(sk_b); random_sk(sk_c); random_sk(sk_d);
    random_sk(sk_e); random_sk(sk_f); random_sk(sk_g); random_sk(sk_h);

    // form the vector of private key
    for (i = 0; i < N; i++) 
      vsk[i] = set_vector(sk_h[i], sk_g[i], sk_f[i], sk_e[i], sk_d[i], sk_c[i], sk_b[i], sk_a[i]);
    
    // intialize the vector of base curve
    for (i = 0; i < HT_NWORDS; i++) {
      vE.y[i] = VSET1(E[0][i]);
      vE.z[i] = VSET1(E[1][i]);
    }

    // keypair generation
    action(&vPK, vsk, &vE);

    gfp_inv_8x1w(t0, vPK.z);
    gfp_mul_8x1w(t1, vPK.y, t0);
    gfp_mont2num_8x1w(t0, t1);

    // We simply swap the position of (pk_a and pk_b), (pk_c and pk_d),
    // (pk_e and pk_f), (pk_g and pk_h) in the AVX512 vectors. 
    // So now it is (pk_g, pk_h, pk_e, pk_f, pk_c, pk_d, pk_a, pk_b)
    //              (sk_h, sk_g, sk_f, sk_e, sk_d, sk_c, sk_b, sk_a)
    for (i = 0; i < HT_NWORDS; i++) { 
      vPK.y[i] = VPERM(vPK.y[i], 0xB1);
      vPK.z[i] = VPERM(vPK.z[i], 0xB1);
    }

    // shared secret computation
    action(&vSS, vsk, &vPK);

    gfp_inv_8x1w(t0, vSS.z);
    gfp_mul_8x1w(t1, vSS.y, t0);
    gfp_mont2num_8x1w(t0, t1);

    get_channel_8x1w(ss_a, t0, 0);
    get_channel_8x1w(ss_b, t0, 1);
    get_channel_8x1w(ss_c, t0, 2);
    get_channel_8x1w(ss_d, t0, 3);
    get_channel_8x1w(ss_e, t0, 4);
    get_channel_8x1w(ss_f, t0, 5);
    get_channel_8x1w(ss_g, t0, 6);
    get_channel_8x1w(ss_h, t0, 7);

    wrong = memcmp(ss_a, ss_b, HT_NWORDS*sizeof(uint32_t));
    wrong |= memcmp(ss_c, ss_d, HT_NWORDS*sizeof(uint32_t));
    wrong |= memcmp(ss_e, ss_f, HT_NWORDS*sizeof(uint32_t));
    wrong |= memcmp(ss_g, ss_h, HT_NWORDS*sizeof(uint32_t));

    if (wrong) {
      printf("No.%d \x1b[31mNOT PASS!\x1b[0m\n", j);
    }
    else 
      printf("No.%d \x1b[32mPASS!\x1b[0m\n", j);
  }

  puts("*******************************************************************");
}

void timing_action()
{
  uint8_t sk_a[N], sk_b[N], sk_c[N], sk_d[N], sk_e[N], sk_f[N], sk_g[N], sk_h[N];
  htpoint vE, vPK;
  __m512i vsk[N];
  int i;

  // randomly generate eight private keys 
  random_sk(sk_a); random_sk(sk_b); random_sk(sk_c); random_sk(sk_d);
  random_sk(sk_e); random_sk(sk_f); random_sk(sk_g); random_sk(sk_h);

  // initialize AVX512 vector of private key
  for (i = 0; i < N; i++) 
    vsk[i] = set_vector(sk_h[i], sk_g[i], sk_f[i], sk_e[i], sk_d[i], sk_c[i], sk_b[i], sk_a[i]);

  for (i = 0; i < HT_NWORDS; i++) {
    vE.y[i] = VSET1(E[0][i]);
    vE.z[i] = VSET1(E[1][i]);
  }

  uint64_t start_cycles, end_cycles, diff_cycles;

  LOAD_CACHE(action(&vPK, vsk, &vE), 100);
  MEASURE_TIME(action(&vPK, vsk, &vE), 1000);
  printf("* GROUP ACTION : %ld cycles\n", diff_cycles);
}

int main() 
{
  test_action();
  // test_multi_actions(1000);
  timing_action();
  
  return 0;
}