#define _XOPEN_SOURCE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <math.h>

#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"
#include "precomputed.h"



GEN torsion_fm;
unsigned long p_length = 256;


GEN norm0(GEN x) {
    return algnorm(global_setup.B, x,0);
}



int test_lideal1() {
    GEN B_O0 = global_setup.O0;
    GEN B = global_setup.B;

    GEN gen = mkcol4s(2,7,3,-5);
    GEN lideal = lideal_create(B, B_O0, gen, NULL);

    //output(lideal);
    output(norm0(gen));
    output(lideal_norm(lideal));


    gel(lideal,3) = NULL;
    GEN z = lideal_generator(lideal);

    output(norm0(z));
    GEN lideal2 = lideal_create(B, B_O0, z, lideal_norm(lideal));

    //output(lideal2);
    gel(lideal2,2) = gen_0;
    output(lideal_norm(lideal2));

    return 0;
}



int test_lideal2() {
    GEN B_O0 = global_setup.O0;
    GEN B = global_setup.B;

    GEN gen = lattice_random(B, B_O0, stoi(10000)); // = mkcol4s(2,7,0,0);

    output(algnorm(B, gen, 0));

    GEN lideal = lideal_create(B, B_O0, gen, NULL);

    printf("LLL\n");
    //output(lideal_lll(lideal));

    //printf("Short\n");
    output(norm0(gel(lideal_lll(lideal),1)));
    output(norm0(gel(lideal_lll(lideal),2)));
    output(norm0(gel(lideal_lll(lideal),3)));
    output(norm0(gel(lideal_lll(lideal),4)));

    //lideal_equiv_prime(lideal);


    GEN I1 = lideal_random_2e(B, B_O0, 256);

    output(lideal_norm(I1));

    printf("A\n");
    GEN I2 = lideal_equiv_prime(I1,NULL);

    output(lideal_norm(I2));


    printf("B\n");
    lideal_isom(I1, I2);

    printf("C\n");
    output(algmultable(B));

    return 0;
}

int test_norm_eq() {
    float accumulated_time_ms = 0.;
    int repetitions = 10000;
    clock_t t;

    //GEN fm = famat_sqr(torsion_fm);

    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        long ctr = 1;
        GEN sol = NULL, target;
        //GEN N = randomi(mpshift(gen_1, 128));
        //GEN p_div_N = gadd(truedivii(global_setup.p, N), gen_1);
        //GEN fm_1;

        t = tic();
        while (!sol) {
            //fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/10) + 2));

            //target = gmul(N,famat_prod(fm_1));
            target = randomi(gmul(global_setup.p, stoi(ctr/20 + 1)));
            sol = norm_equation_special(global_setup.p, target, 0, false);
            ctr++;
        }
        //output(sol);
        accumulated_time_ms += toc(t);

        GEN sum1 = gadd(gsqr(gel(sol,1)), gsqr(gel(sol,2)));
        GEN sum2 = gmul(global_setup.p,gadd(gsqr(gel(sol,3)), gsqr(gel(sol,4))));
        GEN sum = gadd(sum1,sum2);
        if (gcmp(sum,target) != 0) {
            printf("error in norm_equation_special: incorrect result!\n");
            break;
        }

        avma = av;
    }
    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    return 0;
}

int test_famat_rand() {
    float accumulated_time_ms = 0.;
    int repetitions = 100;
    clock_t t;

    GEN fm = famat_sqr(torsion_fm);

    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN target;
        GEN N = randomi(mpshift(gen_1, 128));
        GEN B = gmul(gadd(truedivii(global_setup.p, N), gen_1),stoi(3));
        GEN fm_1;

        t = tic();
        for (int i = 0; i < 100; ++i) {
            //fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/10) + 2));
            fm_1 = famat_random(fm, B,0);
            target = gmul(N,famat_prod(fm_1));
        }
        //output(sol);
        accumulated_time_ms += toc(t);

        avma = av;
    }
    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    return 0;
}

int test_famat() {
    GEN fm = mkmat2(mkcol4s(2,3,5,7), mkcol4s(4,2,1,2));
    GEN fm2;
    output(famat_pop(fm,&fm2));
    output(famat_pop(fm2,&fm2));
    output(famat_pop(fm2,&fm2));
    output(famat_pop(fm2,&fm2));
    output(famat_pop(fm2,&fm2));



    fm = mkmat2(mkcol4s(2,3,5,7), mkcol4s(30,2,3,2));
    output(fm);

    for (int i = 0; i < 20; ++i) {
        fm2 = famat_random(fm,stoi(10),0);
        output(fm2);
    }

    return 0;
}

int test_klpt() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 100;
    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN I = lideal_random_2e(global_setup.B, global_setup.O0, p_length);

        GEN fm = famat_sqr(torsion_fm);

        clock_t t = tic();
        GEN J = klpt_special_smooth(I, fm);
        accumulated_time_ms += toc(t);
        //accumulated_time_ms += toc(t);

        GEN NJ = lideal_norm(J);
        accumulated_bitlength += dbllog2r(itor(NJ,10));
        //printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));

        // check norm

        int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
        if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

        // check isomorphism

        GEN alpha = lideal_isom(I, J); // I*alpha = J
        if (!alpha) { printf("output of klpt is not isomorphic to input\n"); break; }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / repetitions));

    return 0;
}

int test_equiv_nearprime() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 100;
    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN I = lideal_random_2e(global_setup.B, global_setup.O0, p_length);
        GEN fm = famat_sqr(torsion_fm);

        clock_t t = tic();
        GEN J = lideal_equiv_nearprime(I,fm,0);
        // GEN J = lideal_equiv_prime_except(I,NULL,NULL);
        accumulated_time_ms += toc(t);

        //GEN NJ = lideal_norm(J);
        //if (!ispseudoprime(NJ,0)) { printf("output of lideal_equiv_prime_except is not of prime norm\n"); break; }

        // check isomorphism

        GEN alpha = lideal_isom(I, J); // I*alpha = J
        if (!alpha) { printf("output of lideal_equiv_* is not isomorphic to input\n"); break; }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / repetitions));

    return 0;
}

int test_klpt2e() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 100;
    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN I = lideal_random_2e(global_setup.B, global_setup.O0, 64);
        GEN fm = famat_sqr(torsion_fm);

        clock_t t = tic();
        GEN J = klpt_special_smooth_small_2e_input(I, fm);
        accumulated_time_ms += toc(t);
        //accumulated_time_ms += toc(t);

        if (J) {
            GEN NJ = lideal_norm(J);
            accumulated_bitlength += dbllog2r(itor(NJ,10));
            //printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));

            // check norm

            int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
            if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

            // check isomorphism

            GEN alpha = lideal_isom(I, J); // I*alpha = J
            if (!alpha) { printf("output of klpt is not isomorphic to input\n"); }
        }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / repetitions));

    return 0;
}

int test_klpt_general() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 100;
    pari_sp av = avma;
    int win = 0;

    for (int i = 0; i < repetitions; ++i) {
        unsigned int length_NI = 64;

        GEN NI = NULL;

        do {
            NI = randomprime(powiu(gen_2, length_NI));
        } while (Fp_issquare(gen_2,NI));

        GEN alpha = NULL;

        unsigned int margin = 6;
        while (!alpha) {
            alpha = norm_equation_special(global_setup.p, gmul(NI,randomi(powiu(gen_2, p_length-length_NI + margin))), 0, false);
            ++margin;
        }

        GEN I = lideal_create(global_setup.B, global_setup.O0, gtrans(alpha), NI);

        GEN K = lideal_random_2e(global_setup.B, global_setup.O0, p_length);
        K = lideal_equiv_prime(K,NULL);
        clock_t t = tic();
        GEN J = klpt_general_power(I, K, gen_2);
        accumulated_time_ms += toc(t);
        TOC(t,"");


        if (J) {
            ++win;

            GEN NJ = lideal_norm(J);
            accumulated_bitlength += dbllog2r(itor(NJ,10));

            // printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));


            // check norm

            //int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
            //if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

            //output(Z_factor_limit(NJ, 3));

            // check isomorphism

            GEN I1 = lideal_inter(I,K);
            GEN I2 = lideal_inter(I,J);

            GEN alpha = lideal_isom(I1, I2); // I1*alpha = I2
            if (!alpha) { printf("output of klpt is not isomorphic to input! stopped at i=%d \n",i); }
        }

        avma = av;
    }

    printf("klpt general \n");
    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

    return 0;
}

int repres_integer_path_analysis() {
  int repetitions = 50000;
  int global_rep = 1;
  pari_sp av = avma;
  unsigned int length_NI = p_length;
  int size = 1;
  int number = 3* pow(2,size-1);
  GEN pow = gpowgs(gen_2,size);

  printf("size of walk = %d, number of paths = %d \n",size, number);
  printf("number of repetition = %d \n",repetitions);
  GEN NI = NULL;
  GEN list_gen = mklist();
  GEN K1 = lideal_random_2e(global_setup.B, global_setup.O0, size+10);
  K1 = lideal_create(global_setup.B,global_setup.O0,lideal_generator(K1),pow);
  listput(list_gen,lideal_generator(K1),1);
  int list_size =1;
  while (list_size!=number) {
    GEN K2 = lideal_random_2e(global_setup.B, global_setup.O0, size+10);
    K2 = lideal_create(global_setup.B,global_setup.O0,lideal_generator(K2),pow);
    bool test = true;
    for (int i=0;i<list_size;i++){
      GEN K_test = lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],i+1),pow);
      if (lideal_equals(K2,K_test)){
        test =false; break;
      }
    }
    if (test){
      listput(list_gen,lideal_generator(K2),list_size+1);
      // output(list_gen);
      list_size++;
    }
  }
  if (size ==1 ){
    GEN K_spec = lideal_create(global_setup.B,global_setup.O0,mkcol4(gen_1,gen_1,gen_0,gen_0),pow);
    if (lideal_equals(K_spec,lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],1),pow))) {
      printf("the principal ideal of norm 2 is the first \n");
    }
    else if (lideal_equals(K_spec,lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],2),pow))) {
      printf("the principal ideal of norm 2 is the second \n");
    }
    else if (lideal_equals(K_spec,lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],3),pow))) {
      printf("the principal ideal of norm 2 is the third \n");
    }
  }
  int xtot[number];

  for (int i=0;i<number;i++){
    xtot[i] = 0;
  }

  for (int ind = 0;ind < global_rep; ind++){
    int xs[number];

    for (int i=0;i<number;i++){
      xs[i] = 0;
    }




    for (int i = 0; i < repetitions; ++i) {
      pari_sp mid= avma;
      do {
          NI = randomprime(powiu(gen_2, length_NI));
      } while (Fp_issquare(gen_2,NI));

      GEN alpha = NULL;

      unsigned int margin = 20;
      GEN elle = powiu(gen_2, p_length-length_NI + margin);
      // printf("a \n");
      alpha = norm_equation_special_max_order(global_setup.p, gmul(NI,elle), false, true);
      if (alpha) {
        GEN X1;
        alglatcontains(global_setup.B,global_setup.O0,alpha,&X1);
        // output(X1);
        GEN n1 =content(X1);
        // output(n1);
        alpha = mkcol4( gdiv(gel(alpha,1),n1),gdiv(gel(alpha,2),n1),gdiv(gel(alpha,3),n1),gdiv(gel(alpha,4),n1));

        GEN I = lideal_create(global_setup.B, global_setup.O0, alpha, elle);
        GEN I_test = lideal_create(global_setup.B,global_setup.O0,lideal_generator(I),pow);


        for (int j=0;j<number;j++){
          // GEN gen_temp = listpop(list_gen,3);
          GEN J_temp = lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],j+1),pow);
          // listput(list_gen,gen_temp,1);
          if (lideal_equals(J_temp,I_test)){
            xs[j]++;

          }
        }
      }




      avma = mid;
    }
    for (int i=0;i<number;i++){
      xtot[i] += xs[i];
      printf("%d, ",xs[i]);
    }
    printf("\n");
  }
  // printf("total = \n");
  // for (int i=0;i<number;i++){
  //   printf("%d, ",xtot[i]);
  // }
  // printf("\n");



  // printf("x1= %d, x2= %d, x3= %d \n",xs[0],xs[1],xs[2]);
  avma = av;
  // printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
  // printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

  return 0;
}


int sqisign_path_analysis() {
  float accumulated_time_ms = 0., accumulated_bitlength = 0.;
  int repetitions = 50;
  int global_rep = 10;
  pari_sp av = avma;
  int win = 0;
  unsigned int length_NI = 64;
  int size = 3;
  int number = 3* pow(2,size-1);
  GEN pow = gpowgs(gen_2,size);

  printf("size of walk = %d, number of paths = %d \n",size, number);
  printf("number of keys = %d, number of repetition = %d \n",global_rep,repetitions);
  GEN NI = NULL;
  GEN list_gen = mklist();
  GEN K1 = lideal_random_2e(global_setup.B, global_setup.O0, size+10);
  K1 = lideal_create(global_setup.B,global_setup.O0,lideal_generator(K1),pow);
  listput(list_gen,lideal_generator(K1),1);
  int list_size =1;
  while (list_size!=number) {
    GEN K2 = lideal_random_2e(global_setup.B, global_setup.O0, size+10);
    K2 = lideal_create(global_setup.B,global_setup.O0,lideal_generator(K2),pow);
    bool test = true;
    for (int i=0;i<list_size;i++){
      GEN K_test = lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],i+1),pow);
      if (lideal_equals(K2,K_test)){
        test =false; break;
      }
    }
    if (test){
      listput(list_gen,lideal_generator(K2),list_size+1);
      // output(list_gen);
      list_size++;
    }
  }

  if (size ==1 ){
    GEN K_spec = lideal_create(global_setup.B,global_setup.O0,mkcol4(gen_1,gen_1,gen_0,gen_0),pow);
    if (lideal_equals(K_spec,lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],1),pow))) {
      printf("the principal ideal of norm 2 is the first \n");
    }
    else if (lideal_equals(K_spec,lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],2),pow))) {
      printf("the principal ideal of norm 2 is the second \n");
    }
    else if (lideal_equals(K_spec,lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],3),pow))) {
      printf("the principal ideal of norm 2 is the third \n");
    }
  }
  int xtot[number];

  for (int i=0;i<number;i++){
    xtot[i] = 0;
  }

  for (int ind = 0;ind < global_rep; ind++){
    int xs[number];

    for (int i=0;i<number;i++){
      xs[i] = 0;
    }
    do {
        NI = randomprime(powiu(gen_2, length_NI));
    } while (Fp_issquare(gen_2,NI));

    GEN alpha = NULL;

    unsigned int margin = 6;
    while (!alpha) {
        alpha = norm_equation_special(global_setup.p, gmul(NI,randomi(powiu(gen_2, 256-length_NI + margin))), 0, false);
        ++margin;
    }
    GEN I = lideal_create(global_setup.B, global_setup.O0, gtrans(alpha), NI);

    for (int i = 0; i < repetitions; ++i) {
        pari_sp mid= avma;
        GEN K = lideal_random_2e(global_setup.B, global_setup.O0, 130);
        K = lideal_equiv_prime(K,NULL);
        clock_t t = tic();
        GEN J = klpt_general_power(I, K, gen_2);
        accumulated_time_ms += toc(t);
        GEN beta = lideal_generator(J);
        GEN Xb,nb;
        alglatcontains(global_setup.B, global_setup.O0,beta,&Xb);
        nb =content(Xb);
        if (J) {
            ++win;

            GEN NJ = lideal_norm(J);
            accumulated_bitlength += dbllog2r(itor(NJ,10));

            GEN J_test = lideal_create(global_setup.B,global_setup.O0,lideal_generator(J),pow);


            for (int j=0;j<number;j++){
              // GEN gen_temp = listpop(list_gen,3);
              GEN J_temp = lideal_create(global_setup.B,global_setup.O0,gel(list_gen[2],j+1),pow);
              // listput(list_gen,gen_temp,1);
              if (lideal_equals(J_temp,J_test)){
                xs[j]++;
                break;
              }
            }

            GEN I1 = lideal_inter(I,K);
            GEN I2 = lideal_inter(I,J);

            GEN alpha = lideal_isom(I1, I2); // I1*alpha = I2
            if (!alpha) { printf("output of klpt is not isomorphic to input! stopped at i=%d \n",i); }
        }


        avma = mid;
    }
    int checksum =0;
    for (int i=0;i<number;i++){
      xtot[i] += xs[i];
      checksum+=xs[i];
      printf("%d, ",xs[i]);
    }
    printf("\n");
    if (checksum!= repetitions) {
        printf("the total is off !!\n");
    }
  }
  printf("total = \n");
  for (int i=0;i<number;i++){
    printf("%d, ",xtot[i]);
  }
  printf("\n");



  // printf("x1= %d, x2= %d, x3= %d \n",xs[0],xs[1],xs[2]);
  avma = av;
  // printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
  // printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

  return 0;
}



int test_eichler_order_norm_equation() {
  float accumulated_time_ms = 0., accumulated_bitlength = 0.;
  int repetitions = 100;
  pari_sp av = avma;
  int win = 0;
  // printf("wow \n");

  for (int i = 0; i < repetitions; ++i) {
      GEN l=gen_2;
      unsigned int length_NI = 128;

      GEN NI = NULL;

      do {
          NI = randomprime(powiu(gen_2, length_NI));
      } while (Fp_issquare(gen_2,NI));

      GEN alpha = NULL;

      unsigned int margin = 6;
      while (!alpha) {
          alpha = norm_equation_special(global_setup.p, gmul(NI,randomi(powiu(gen_2, 256-length_NI + margin))), 0, false);
          ++margin;
      }
      // printf("bloopers \n");
      GEN test_alg = mkcol4(gdiv(gen_1,gen_2),gen_1,gdiv(gen_1,gen_2),gen_1);
      // output(test_alg);
      // output(global_setup.O0);
      // if (!alglatcontains(global_setup.B,global_setup.O0,test_alg,NULL)) {
      //   printf("not in \n");
      // }
      // GEN I = lideal_create(global_s  printf("eichler norm equation \n");etup.B, global_setup.O0, gtrans(alpha), gdivexact(algnorm(global_setup.B,alpha,0),NI));
      GEN I = lideal_random_2e(global_setup.B, global_setup.O0,256);
      GEN fm = famat_sqr(torsion_fm);


      // output(gel())
      // output(gel(fm,1));
      // gel(gel(fm,2),1) = gen_0; // remove power of 3
      GEN L = lideal_random_2e(global_setup.B,lideal_right_order(I),65);
      eichler_package e;
      e.beta=NULL; e.J = NULL; e.delta=NULL; e.I = NULL; e.L=NULL;



      clock_t t = tic();
      eichler_norm_equation_special_smooth(&e,I, L, fm, l);
      float tim = toc(t);
      // printf("%f ms \n",tim);
      accumulated_time_ms += tim;
      if (e.J) {
          ++win;
          GEN No=algnorm(global_setup.B,e.beta,0);
          GEN NJ = lideal_norm(e.J);
          accumulated_bitlength += dbllog2r(itor(No,10));



          // check norm

          //int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
          //if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

          //output(Z_factor_limit(NJ, 3));

          //check if the solution is good. Start by verifying that delta is in the correct eichler order
          GEN beta=e.beta;
          // GEN beta_bar= alg_conj(global_setup.B,beta);


          GEN Tr=algtrace(global_setup.B,beta,0);
          // GEN lpow=lideal_norm(L);
          Tr=Fp_mul(Tr,Fp_inv(gen_2,NJ),NJ);
          GEN t=Fp_sub( Fp_mul(Tr,Tr,NJ),No,NJ );
          GEN sqN=Fp_sqrt(t,NJ);
          if (!sqN) {printf("no solution to the sqrt mod NJ\n");}
          Tr= gneg(Tr);
          GEN x1N= gadd(Tr,sqN);
          GEN x2N=gsub(Tr,sqN);
          x1N=mkcol4(x1N,gen_0,gen_0,gen_0);
          x2N=mkcol4(x2N,gen_0,gen_0,gen_0);
          // output(x1N);
          // x1N=gmul(x1N,1);
          // a=[x1N,0,0,0]~;
          // x1N=algtomatrix(global_setup.B,a,1);
          // output(x1N);
          // x1N=alglathnf(global_setup.B,x1N,gen_0);
          //

          // x2N=algtomatrix(global_setup.B,mkcol4(x2N,0,0,0),1);
          //
          // x2N=alglathnf(global_setup.B,x2N,gen_0);
          if (! (
            alglatcontains(global_setup.B,gel(e.J,1), algadd(global_setup.B,beta,x1N ),NULL )
            || alglatcontains(global_setup.B,gel(e.J,1), algadd(global_setup.B,beta,x2N ),NULL )
          ) ) { printf(" not in the eichler order \n "); }

          if (!alglatcontains(lideal_algebra(e.J),gel(e.J,1),alg_conj(lideal_algebra(e.J),e.delta),NULL )) {
            printf("delta bar is not in the ideal J\n");
          }

          //second part of the verification with find_coeff
          GEN A = global_setup.B;
          GEN K = lideal_random_2e(global_setup.B,lideal_order(L),33);
          // printf("before coeff 2\n");
          GEN coeffs = find_coeff(K,L,beta,e.delta, lideal_right_order(e.J));
          // printf("after coeff \n");
          if (coeffs) {
            GEN C0 = gel(coeffs,1);
            GEN D0 = gel(coeffs,2);
            // output(C0);
            // output(D0);
            GEN C0_alg = mkcol4(C0,gen_0,gen_0,gen_0);
            GEN D0_alg = mkcol4(D0,gen_0,gen_0,gen_0);

            GEN new_beta = algmul(A,beta,e.delta);
            new_beta = algmul(A,alg_conj(A,e.delta),new_beta);
            GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,e.delta,0)),gen_0,gen_0,gen_0);
            new_beta = algmul(A,new_beta,overnorm);


            GEN theta = algadd(A,algmul(A,D0_alg,new_beta),C0_alg);
            GEN theta_bar = alg_conj(A,theta);
            GEN alpha = lideal_generator(L);


            if (!alglatcontains(A,lideal_lattice(K), algmul(A,alpha,theta_bar),NULL )) {
              printf("it is not in \n");
            }
            GEN K_test= lideal_create(A,lideal_order(K),algmul(A,alpha,theta_bar),lideal_norm(K));




            int check2 = lideal_equals(K_test,K);
            if (!check2) {

              printf("didn't work \n");
              output(C0);
              output(D0);
              // output(vali(gen_1))
              // printf("%d \n",Z_pval(lideal_norm(K),gen_2));
              // output(Z_pval(lideal_norm(K),gen_2));
              printf("%ld \n",vali( algnorm(A,theta_bar,0) ));
              printf("%ld \n",vali(lideal_norm(K_test)));
              printf("%ld \n",vali(algnorm(A,algmul(A,alpha,theta_bar),0)));
            }
          }
          else {
            printf("no solution was found for coeffs \n");
          }
          // GEN coeffs2 = find_coeff2(K,L,beta,e.delta, lideal_right_order(e.J));
          //
          // GEN I1 = lideal_inter(I,K);
          // GEN I2 = lideal_inter(I,J);
          //
          // GEN alpha = lideal_isom(I1, I2); // I1*alpha = I2
          // if (!alpha) { printf("output of klpt is not isomorphic to input!\n"); break; }
      }
      avma = av;
  }
  printf("eichler norm equation \n");
  printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
  printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

  return 0;
}



int test_eichler_order_norm_equation_small() {
  float accumulated_time_ms = 0., accumulated_bitlength = 0.;
  int repetitions = 30;
  pari_sp av = avma;
  int win = 0;


  for (int i = 0; i < repetitions; ++i) {
      GEN l=gen_2;
      unsigned int length_NI = 128;

      GEN NI = NULL;

      do {
          NI = randomprime(powiu(gen_2, length_NI));
      } while (Fp_issquare(gen_2,NI));


      GEN I = lideal_random_2e(global_setup.B, global_setup.O0,256);
      GEN fm = famat_sqr(torsion_fm);
      // output(gel())
      // output(gel(fm,1));
      // gel(gel(fm,2),1) = gen_0; // remove power of 3
      GEN L = lideal_random_2e(global_setup.B,lideal_right_order(I),33);
      eichler_package e;
      e.beta=NULL; e.J = NULL; e.delta=NULL; e.I = NULL; e.L=NULL;
      GEN remove_fact = mkmat2(mkcol4s(6983,4019,3691,4283),mkcol4s(2,2,2,2));
      fm = famat_div(fm,remove_fact);
      remove_fact = mkmat2(mkcols(2713),mkcols(2));
      fm = famat_div(fm,remove_fact);



      clock_t t = tic();
      eichler_norm_equation_special_smooth_small_curve(&e,I, L, fm, l);
      float tim = toc(t);
      // printf("%f ms \n",tim);
      accumulated_time_ms += tim;
      if (e.J) {
          ++win;
          GEN No=algnorm(global_setup.B,e.beta,0);
          GEN NJ = lideal_norm(e.J);

          accumulated_bitlength += dbllog2r(itor(No,10));



          // check norm

          //int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
          //if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

          //output(Z_factor_limit(NJ, 3));

          //check if the solution is good. Start by verifying that delta is in the correct eichler order
          GEN beta=e.beta;
          // GEN beta_bar= alg_conj(global_setup.B,beta);



          if (!alglatcontains(global_setup.B,lideal_right_order(e.J),e.beta,NULL )) {
            printf("beta is not in the order OR(J)\n");
          }


          GEN Tr=algtrace(global_setup.B,beta,0);
          // GEN lpow=lideal_norm(L);
          Tr=Fp_mul(Tr,Fp_inv(gen_2,NJ),NJ);
          GEN t=Fp_sub( Fp_mul(Tr,Tr,NJ),No,NJ );
          GEN sqN=Fp_sqrt(t,NJ);
          if (!sqN) {printf("no solution to the sqrt mod NJ\n");}
          Tr= gneg(Tr);
          GEN x1N= gadd(Tr,sqN);
          GEN x2N=gsub(Tr,sqN);
          x1N=mkcol4(x1N,gen_0,gen_0,gen_0);
          x2N=mkcol4(x2N,gen_0,gen_0,gen_0);
          // output(x1N);
          // x1N=gmul(x1N,1);
          // a=[x1N,0,0,0]~;
          // x1N=algtomatrix(global_setup.B,a,1);
          // output(x1N);
          // x1N=alglathnf(global_setup.B,x1N,gen_0);
          //

          // x2N=algtomatrix(global_setup.B,mkcol4(x2N,0,0,0),1);
          //
          // x2N=alglathnf(global_setup.B,x2N,gen_0);
          if (! (
            alglatcontains(global_setup.B,gel(e.J,1), algadd(global_setup.B,beta,x1N ),NULL )
            || alglatcontains(global_setup.B,gel(e.J,1), algadd(global_setup.B,beta,x2N ),NULL )
          ) ) { printf(" not in the eichler order \n "); }

          if (!alglatcontains(lideal_algebra(e.J),gel(e.J,1),alg_conj(lideal_algebra(e.J),e.delta),NULL )) {
            printf("delta bar is not in the ideal J\n");
          }

          //second part of the verification with find_coeff
          GEN A = global_setup.B;
          GEN K = lideal_random_2e(global_setup.B,lideal_order(L),33);
          // printf("before coeff 2\n");
          GEN coeffs = find_coeff(K,L,beta,e.delta, lideal_right_order(e.J));
          // printf("after coeff \n");
          if (coeffs) {
            GEN C0 = gel(coeffs,1);
            GEN D0 = gel(coeffs,2);
            // output(C0);
            // output(D0);
            GEN C0_alg = mkcol4(C0,gen_0,gen_0,gen_0);
            GEN D0_alg = mkcol4(D0,gen_0,gen_0,gen_0);

            GEN new_beta = algmul(A,beta,e.delta);
            new_beta = algmul(A,alg_conj(A,e.delta),new_beta);
            GEN overnorm = mkcol4(gdiv(gen_1,algnorm(A,e.delta,0)),gen_0,gen_0,gen_0);
            new_beta = algmul(A,new_beta,overnorm);


            GEN theta = algadd(A,algmul(A,D0_alg,new_beta),C0_alg);
            GEN theta_bar = alg_conj(A,theta);
            GEN alpha = lideal_generator(L);


            if (!alglatcontains(A,lideal_lattice(K), algmul(A,alpha,theta_bar),NULL )) {
              printf("it is not in \n");
            }
            GEN K_test= lideal_create(A,lideal_order(K),algmul(A,alpha,theta_bar),lideal_norm(K));




            int check2 = lideal_equals(K_test,K);
            if (!check2) {

              printf("didn't work \n");
              output(C0);
              output(D0);
              // output(vali(gen_1))
              // printf("%d \n",Z_pval(lideal_norm(K),gen_2));
              // output(Z_pval(lideal_norm(K),gen_2));
              printf("%ld \n",vali( algnorm(A,theta_bar,0) ));
              printf("%ld \n",vali(lideal_norm(K_test)));
              printf("%ld \n",vali(algnorm(A,algmul(A,alpha,theta_bar),0)));
            }
          }
          else {
            printf("no solution was found for coeffs \n");
          }
          // GEN coeffs2 = find_coeff2(K,L,beta,e.delta, lideal_right_order(e.J));
          //
          // GEN I1 = lideal_inter(I,K);
          // GEN I2 = lideal_inter(I,J);
          //
          // GEN alpha = lideal_isom(I1, I2); // I1*alpha = I2
          // if (!alpha) { printf("output of klpt is not isomorphic to input!\n"); break; }
      }
      avma = av;
  }
  printf("eichler norm equation with next extremal order \n");
  printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
  printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

  return 0;
}

// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(80000000, 1<<18);
    init_precomputations();

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }
    torsion_fm = famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact);

    // test_famat();
    // test_klpt();
    test_klpt_general();
    // test_eichler_order_norm_equation() ;
    // test_eichler_order_norm_equation_small() ;
    // sqisign_path_analysis();
    // repres_integer_path_analysis();

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);
}
