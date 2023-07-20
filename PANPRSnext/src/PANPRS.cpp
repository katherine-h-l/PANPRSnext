#include <RcppArmadillo.h>

#include <stdio.h>
#include <stdlib.h>

#include "PANPRS.hpp"

Rcpp::List gsfPEN_cpp(
    arma::Mat<double> summary_betas,
    arma::Col<int> ld_J,
    arma::Mat<int> index_matrix,
    arma::Col<int> index_J,
    arma::Col<double> ld_vec,
    arma::Mat<double> SD_vec,
    arma::Mat<double> tuning_matrix,
    arma::Col<double> lambda0_vec,
    arma::Mat<double> z_matrix,
    arma::Col<double> lambda_vec,
    arma::Mat<int> func_lambda,
    arma::Col<int> Ifunc_SNP,
    arma::Col<int> dims,
    arma::Col<double> params)
{
  int num_SNP = dims(0);

  int P = dims(1);
  int Q = dims(2);

  int nrow_index_matrix = dims(3);
  int ncol_index_matrix = dims(4);

  int nrow_z_matrix = dims(5);
  int ncol_z_matrix = dims(6);

  int nrow_func_lambda = dims(7);
  int ncol_func_lambda = dims(8);

  int nrow_tuning_matrix = dims(9);
  int ncol_tuning_matrix = dims(10);

  int nrow_all_tuning_matrix = dims(11);
  int ncol_all_tuning_matrix = dims(12);

  int nrow_beta_matrix = dims(13);
  int ncol_beta_matrix = dims(14);

  double upper_val = params(0);
  int num_iter = params(1);
  int breaking = params(2);
  int z_scale = params(3);
  int df_max = params(4);
  int leng_p_threshold = params(5);
  int num_indices = params(6);

  double epsilon = 0.0001;

  arma::Col<int> num_iter_vec(nrow_all_tuning_matrix, arma::fill::zeros);
  arma::Col<double> temp_lambda_vec(ncol_func_lambda, arma::fill::zeros);
  arma::Col<double> sum_betas(P, arma::fill::zeros);

  arma::Mat<double> all_tuning_matrix(nrow_all_tuning_matrix, ncol_all_tuning_matrix, arma::fill::zeros);
  // arma::Mat<double> beta_matrix(nrow_beta_matrix, ncol_beta_matrix, arma::fill::zeros);
  arma::SpMat<double> beta_matrix(nrow_beta_matrix, ncol_beta_matrix);
  arma::Mat<double> joint_b_matrix(P, Q, arma::fill::zeros);
  arma::Mat<double> temp_b_matrix(P, Q, arma::fill::zeros);
  arma::Mat<int> skip(P, Q, arma::fill::zeros);

  GetRNGstate();

  // len_p_Thresold = 4
  // By default, p.Threshold = seq(0.5, 10^-4, length.out=4)
  int max_tuning_index = nrow_all_tuning_matrix * nrow_func_lambda * nrow_tuning_matrix;
  int tuning_index = -1;
  for (int threshold_index = 0; threshold_index < leng_p_threshold; threshold_index++)
  {
    for (int tun_idx_1 = 0; tun_idx_1 < nrow_func_lambda; tun_idx_1++)
    {
      for (int tun_idx_2 = 0; tun_idx_2 < nrow_tuning_matrix; tun_idx_2++)
      {
        tuning_index++;

        for (int func_index = 0; func_index < ncol_func_lambda; func_index++)
        {
          temp_lambda_vec(func_index) = lambda_vec(func_lambda(tun_idx_1, func_index));
          all_tuning_matrix(tuning_index, func_index + 1) = temp_lambda_vec(func_index);
        }

        // lambda0vec = abs(-qnorm(p.Threshold/2)) -> 0.67 0.96 1.38 3.89
        double lambda0 = lambda0_vec(threshold_index);
        all_tuning_matrix(tuning_index, 0) = lambda0;

        double lambda2 = tuning_matrix(tun_idx_2, 2);
        all_tuning_matrix(tuning_index, ncol_all_tuning_matrix - 2) = lambda2;

        double tau2 = tuning_matrix(tun_idx_2, 3);
        all_tuning_matrix(tuning_index, ncol_all_tuning_matrix - 1) = tau2;

        bool converges = true;
        for (int n = 1; n <= num_iter; n++)
        {
          // This is non zero when there are Z-scores (summary_betas) that are not in the linkage data (LdJ)
          // Realistically, this is always 0 given the R interfaces only selects those that are present in both
          if (num_indices != 0)
          {
            // num_indices = 3, 1: SNP_A, 2: SNP_B, 3: R
            for (int i = 0; i < num_indices; i++)
            {
              int j = index_J(i);
              double lambda1 = lambda0;

              if (Ifunc_SNP(j) == 1)
              {
                for (int func_index = 0; func_index < ncol_func_lambda; func_index++)
                {
                  lambda1 = lambda1 + z_matrix(j, func_index) * temp_lambda_vec(func_index);
                }
              }

              for (int q = 0; q < Q; q++)
              {
                double bj_bar = summary_betas(j, q);
                if (bj_bar != 0.0)
                {
                  double threshold;
                  if (Q == 1)
                    threshold = lambda1;
                  else
                    threshold = lambda1 + lambda2 / (sum_betas(j) + tau2);

                  // SDvec is 1/sqrt(sample size of curr. GWA)
                  if (z_scale == 1)
                    threshold *= SD_vec(j, q);

                  if (bj_bar > threshold)
                    joint_b_matrix(j, q) = bj_bar - threshold;
                  else if (bj_bar < -threshold)
                    joint_b_matrix(j, q) = bj_bar + threshold;
                  else
                    joint_b_matrix(j, q) = 0.0;
                }

                if (summary_betas(j, q) * joint_b_matrix(j, q) < 0)
                {
                  Rprintf("summary_betas[j]=%d\n", j);
                  Rprintf("summary_betas[q]=%d\n", q);
                  Rprintf("summary_betas[j][q]=%e\n", summary_betas(j, q));
                  Rprintf("joint_b_matrix[j][q]=%e\n", joint_b_matrix(j, q));
                  perror("sign inverse");
                }
              }
            }
          }

          // Everything up to this point is only done on the first iteration, of course
          // once for each set of tuning params (standard and functional)
          // Inserting a “#pragma omp parallel for” right here (with some other boilerplate) should be
          // All that is necessary to parallelize it (to start) (will test this)
          // Since there are lots of nested for loops, it may be possible to parallelize it further with
          // “collapse(n)” where n is the number of nested loops to parallelize
          for (int p = 0; p < num_SNP; p++)
          {
            int j = index_matrix(p, 0);
            double lambda1 = lambda0;

            if (Ifunc_SNP(j) == 1)
            {
              for (int func_index = 0; func_index < ncol_func_lambda; func_index++)
              {
                lambda1 += z_matrix(j, func_index) * temp_lambda_vec(func_index);
              }
            }

            for (int q = 0; q < Q; q++)
            {
              if (skip(j, q) == 0)
              {
                if (summary_betas(j, q) != 0.0)
                {
                  double bj_bar = summary_betas(j, q);
                  int i = index_matrix(p, 1);

                  for (int i = index_matrix(p, 1); i < index_matrix(p, 2) + 1; i++)
                  {
                    bj_bar -= ld_vec(i) * joint_b_matrix(ld_J(i), q);
                  }
                  // if (tuning_index == 18 && n == 1) {
                  //   printf("q: %d, i: %d, ld_vec(i): %f\n", q, i, ld_vec(i), joint_b_matrix(ld_J(i), q));
                  //   printf("ld_J(i): %d, joint_b_matrix(ld_J(i), q): %f\n", ld_J(i), joint_b_matrix(ld_J(i), q));
                  //   printf("bj_bar: %f\n", bj_bar);
                  // }

                  // if (tuning_index == 18 && n == 1) printf("lambda1: %f\n", lambda1);
                  double threshold;
                  if (Q == 1)
                    threshold = lambda1;
                  else
                    threshold = lambda1 + lambda2 / (sum_betas(j) + tau2);

                  // if (tuning_index == 18 && n == 1) {
                  //   printf("threshold: %f\n", threshold);
                  //   printf("sum_betas(j): %f\n", sum_betas(j));
                  // }

                  if (fabs(bj_bar) > upper_val)
                  {
                    if (breaking == 1)
                    {
                      num_iter_vec(tuning_index) = -1;
                      break;
                    }
                    else
                    {
                      bj_bar = 0.0;
                      skip(j, q) = 1;
                    }
                  }

                  if (z_scale == 1)
                    threshold *= SD_vec(j, q);

                  if (bj_bar > threshold)
                    joint_b_matrix(j, q) = bj_bar - threshold;
                  else if (bj_bar < -threshold)
                    joint_b_matrix(j, q) = bj_bar + threshold;
                  else
                    joint_b_matrix(j, q) = 0.0;
                }
                else
                {
                  joint_b_matrix(j, q) = 0.0;
                }
              }
              else
              {
                joint_b_matrix(j, q) = 0.0;
              }
            }
          }

          // Actual iteration that checks convergence
          bool found = false;
          int df_q = 0;
          for (int p = 0; p < P; p++)
          {
            for (int q = 0; q < Q; q++)
            {
              if (fabs(joint_b_matrix(p, q)) > upper_val)
                skip(p, q) = 1;

              if (q == 0)
              {
                if (fabs(joint_b_matrix(p, q)) != 0)
                  df_q++;

                if (df_q > df_max)
                {
                  converges = false;
                  break;
                }
              }

              if (fabs(temp_b_matrix(p, q) - joint_b_matrix(p, q)) > epsilon)
              {
                found = true;
                break;
              }
            }
          }

          if (!converges)
          {
            num_iter_vec(tuning_index) = -2;
            break;
          }

          if (!found)
          {
            int beta_idx = 0;
            for (int q = 0; q < Q; q++)
            {
              for (int p = 0; p < P; p++)
              {
                if (joint_b_matrix(p, q) != 0.0)
                  beta_matrix(tuning_index, beta_idx) = joint_b_matrix(p, q);
                beta_idx++;
              }
            }

            num_iter_vec(tuning_index) = n;
            break;
          }

          for (int p = 0; p < num_SNP; p++)
          {
            int j = index_matrix(p, 0);
            sum_betas(j) = 0.0;

            for (int q = 0; q < Q; q++)
            {
              temp_b_matrix(j, q) = joint_b_matrix(j, q);
              sum_betas(j) += fabs(joint_b_matrix(j, q));
            }
          }

          if (num_indices != 0)
          {
            for (int i = 0; i < num_indices; i++)
            {
              int j = index_J(i);
              sum_betas(j) = 0.0;
              for (int q = 0; q < Q; q++)
              {
                temp_b_matrix(j, q) = joint_b_matrix(j, q);
                sum_betas(j) += fabs(joint_b_matrix(j, q));
              }
            }
          }

          num_iter_vec(tuning_index) = n;
        }
      }
    }
  }
  PutRNGstate();

  return Rcpp::List::create(
      Rcpp::Named("beta_matrix") = beta_matrix,
      Rcpp::Named("num_iter_vec") = num_iter_vec,
      Rcpp::Named("all_tuning_matrix") = all_tuning_matrix);
}
