data {
    int<lower=1> row_n;
    int<lower=1> col_n;
    int<lower=1> colcat_n;
    int<lower=0> total[row_n * col_n];
    int<lower=1,upper=row_n> row_index[row_n * col_n];
    int<lower=1,upper=col_n> col_index[row_n * col_n];
    int<lower=1,upper=colcat_n> colcat_index[col_n];
    int<lower=0> row_count[row_n * col_n];
    int<lower=0> col_count[row_n * col_n];
    int<lower=0> cell_count[row_n * col_n];
}
transformed data {
    int<lower=1> n = row_n * col_n;
    int<lower=0> cell_count_compl[n];
    int<lower=0> row_count_compl[n];
    for(i in 1:n){
      cell_count_compl[i] = col_count[i] - cell_count[i];
      row_count_compl[i] = total[i] - row_count[i];
    }
}
parameters {
    //col params
    real col_mean;
    real<lower=0> col_sd;
    vector[col_n] raw_col_logodds;
    real<lower=0> cell_sd;
    vector[n] raw_cell_logodds;

    //biases in deviations terms
    vector[col_n] raw_col_bias;
    vector[n] raw_cell_bias;
    real<lower=0> col_bias_sd;
    real<lower=0> cell_bias_sd;
}
transformed parameters {
    //recenter params
    vector[col_n] col_logodds = raw_col_logodds * col_sd;
    vector[n] cell_logodds = raw_cell_logodds * cell_sd + col_logodds[col_index];

    //incorporate bias
    vector[col_n] col_bias = raw_col_bias * col_bias_sd;
    vector[n] cell_bias = raw_cell_bias * cell_bias_sd;
    vector[n] cell_logodds_focal = cell_logodds +
              (col_bias[col_index] + cell_bias) / 2;
    vector[n] cell_logodds_compl = cell_logodds -
              (col_bias[col_index] + cell_bias) / 2;
}
model {
    //priors and hyperpriors

    //marginal params
    col_mean ~ normal(0,2);
    col_sd ~ std_normal();
    raw_col_logodds ~ std_normal();
    raw_cell_logodds ~ std_normal();
    cell_sd ~ std_normal();

    //bias params
    raw_col_bias ~ std_normal();
    col_bias_sd ~ std_normal();

    raw_cell_bias ~ std_normal();
    cell_bias_sd ~ std_normal();

    //likelihood
    cell_count ~ binomial_logit(row_count, cell_logodds_focal);
    cell_count_compl ~ binomial_logit(row_count_compl, cell_logodds_compl);

}
generated quantities {
    vector[n] cell_total_prob_bias = inv_logit(cell_logodds_focal) - inv_logit(cell_logodds_compl);
}